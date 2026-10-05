"""Evaluate EMERGE rules with expressions and counts, rather than graphs.

There are three distinct uses of the same grammar: KO combinations inside
reactions, numbered reactions inside pathways, and strict signature requirements.
Keeping those scoring policies separate preserves AnnoGuild's original results.
"""

import csv
from dataclasses import dataclass
from itertools import product
from pathlib import Path
import re


DATA_DIR = Path(__file__).parent / "data"


@dataclass(frozen=True)
class Expression:
    """A token, an AND (+), or an OR (,); children retain definition order."""

    kind: str
    token: str = ""
    children: tuple = ()


def parse_expression(definition: str) -> Expression:
    """Parse parentheses, '+' and ','; AND binds more tightly than OR.

    Example: '(K00001+K00002),K00003' is OR(AND(K00001,K00002),K00003).
    No Python evaluation or graph library is involved.
    """
    tokens = re.findall(r"[A-Za-z0-9_:.-]+|[()+,]|\S", definition)
    position = 0

    def parse_group(separator, child_parser, kind):
        nonlocal position
        children = [child_parser()]
        while position < len(tokens) and tokens[position] == separator:
            position += 1
            children.append(child_parser())
        return children[0] if len(children) == 1 else Expression(kind, children=tuple(children))

    def parse_atom():
        nonlocal position
        if position == len(tokens):
            raise ValueError(f"Missing token in rule: {definition!r}")
        token = tokens[position]
        position += 1
        if token == "(":
            expression = parse_or()
            if position == len(tokens) or tokens[position] != ")":
                raise ValueError(f"Unbalanced parentheses: {definition!r}")
            position += 1
            return expression
        if not re.fullmatch(r"[A-Za-z0-9_:.-]+", token):
            raise ValueError(f"Unexpected token {token!r} in rule: {definition!r}")
        return Expression("token", token=token)

    def parse_and():
        return parse_group("+", parse_atom, "and")

    def parse_or():
        return parse_group(",", parse_and, "or")

    expression = parse_or()
    if position != len(tokens):
        raise ValueError(f"Unexpected trailing tokens in rule: {definition!r}")
    return expression


def expression_tokens(expression: Expression) -> set[str]:
    if expression.kind == "token":
        return {expression.token}
    return set().union(*(expression_tokens(child) for child in expression.children))


def signature_present(expression: Expression, features: set[str]) -> bool:
    """Signatures use strict Boolean logic: every AND member is required."""
    if expression.kind == "token":
        return expression.token in features
    calls = [signature_present(child, features) for child in expression.children]
    return all(calls) if expression.kind == "and" else any(calls)


def reaction_options(expression: Expression, features: set[str]) -> list[tuple[str, ...]]:
    """Return alternative KO requirements eligible for legacy coverage scoring.

    Missing AND subunits remain in the denominator. Entirely absent OR alternatives
    are skipped. A reaction is partially present if *any* of its KOs is present;
    strict AND is used only for signatures, not this legacy presence check.
    """
    if expression.kind == "token":
        return [(expression.token,)]
    if expression.kind == "or":
        return [option for child in expression.children
                if expression_tokens(child) & features
                for option in reaction_options(child, features)]
    return [sum(combination, ()) for combination in product(
        *(reaction_options(child, features) for child in expression.children)
    )]


def pathway_routes(expression: Expression) -> list[tuple[str, ...]]:
    """Expand numbered reaction alternatives into explicit pathway routes."""
    if expression.kind == "token":
        return [(expression.token,)]
    children = [pathway_routes(child) for child in expression.children]
    if expression.kind == "or":
        return [route for alternatives in children for route in alternatives]
    return [sum(combination, ()) for combination in product(*children)]


@dataclass(frozen=True)
class PathwayRule:
    pathway_id: str
    reactions: dict[str, Expression]
    routes: tuple
    signature: Expression | None


def read_tsv(path) -> tuple[list[str], list[dict[str, str]]]:
    """Read TSVs without numeric/Boolean coercion or a pandas dependency."""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"Missing or duplicate column names in {path}")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError(f"A row has the wrong number of TSV fields in {path}")
        return reader.fieldnames, rows


def load_rules(reactions_file=None, pathways_file=None) -> list[PathwayRule]:
    """Load explicit reaction IDs and validate all references before evaluation."""
    _, rows = read_tsv(reactions_file or DATA_DIR / "reactions.tsv")
    reactions = {}
    for row in rows:
        key = (row["pathway"], row["reaction_id"])
        if key in reactions:
            raise ValueError(f"Duplicate reaction: {key}")
        expression = parse_expression(row["definition"])
        if any(not re.fullmatch(r"K\d{5,6}", token)
               for token in expression_tokens(expression)):
            raise ValueError(f"Invalid KO/custom feature in reaction {key}")
        reactions[key] = expression

    _, rows = read_tsv(pathways_file or DATA_DIR / "pathways.tsv")
    rules = []
    seen = set()
    for row in rows:
        pathway_id = row["pathway"] + "-" + row["subpathway"]
        if pathway_id in seen:
            raise ValueError(f"Duplicate pathway: {pathway_id}")
        seen.add(pathway_id)
        routes = pathway_routes(parse_expression(row["reaction"]))
        ids = {reaction for route in routes for reaction in route}
        missing = ids - {key[1] for key in reactions if key[0] == row["pathway"]}
        if missing:
            raise ValueError(f"Unknown reactions in {pathway_id}: {sorted(missing)}")
        signature = row["signature_definition"].strip()
        signature_expression = parse_expression(signature) if signature else None
        if signature_expression and any(not re.fullmatch(r"K\d{5,6}", token)
                                        for token in expression_tokens(signature_expression)):
            raise ValueError(f"Invalid signature KO/custom feature in {pathway_id}")
        rules.append(PathwayRule(
            pathway_id,
            {reaction: reactions[row["pathway"], reaction] for reaction in ids},
            tuple(routes), signature_expression,
        ))
    return rules


def evaluate_pathway(rule: PathwayRule, features: set[str],
                     reaction_threshold=0.7, ko_threshold=0.6) -> dict:
    """Score one genome/pathway using AnnoGuild's two independent maxima.

    A completely absent reaction contributes zero to reaction coverage, but its
    KOs are excluded from KO coverage. Repeated KOs in different reactions count
    separately. Keep these unusual legacy choices explicit rather than silently
    changing the scientific method during a refactor.
    """
    reaction_coverage = 0.0
    best_ko_coverage = 0.0
    best_required = ()
    best_route = ()
    for route in rule.routes:
        observed = [bool(expression_tokens(rule.reactions[r]) & features) for r in route]
        reaction_coverage = max(reaction_coverage, sum(observed) / len(route))
        # Dynamic programming stores only distinct (found,total) counts, avoiding
        # enumerating every combination of enzyme alternatives across reactions.
        counts = {(0, 0): ()}
        for reaction, present in zip(route, observed):
            options = reaction_options(rule.reactions[reaction], features) if present else [()]
            updated = {}
            for (found, total), required in counts.items():
                for option in options:
                    key = (found + sum(ko in features for ko in option), total + len(option))
                    updated.setdefault(key, required + option)
            counts = updated
        for (found, total), required in counts.items():
            coverage = found / total if total else 0.0
            if coverage > best_ko_coverage:
                best_ko_coverage, best_required, best_route = coverage, required, route

    signature_ok = rule.signature is None or signature_present(rule.signature, features)
    required = set(best_required)
    # Include all alternatives for entirely missing reactions in the selected
    # route. These are explanatory candidates, not additional KO denominators.
    for reaction in best_route:
        tokens = expression_tokens(rule.reactions[reaction])
        if not tokens & features:
            required.update(tokens)
    if not best_route:
        required = set().union(*(expression_tokens(e) for e in rule.reactions.values()))
    all_required = set().union(*(expression_tokens(e) for e in rule.reactions.values()))
    signature_required = expression_tokens(rule.signature) if rule.signature else set()
    called = (signature_ok and reaction_coverage > 0
              and reaction_coverage >= reaction_threshold and best_ko_coverage >= ko_threshold)
    # Additional reports never change the Boolean product or specialization schema.
    return {
        "pathway": rule.pathway_id,
        "called": called,
        "reaction_coverage": reaction_coverage,
        "ko_coverage": best_ko_coverage,
        "signature_present": signature_ok,
        "ko_coverage_route": "+".join(best_route),
        "supporting_features": ";".join(sorted(all_required & features)),
        "missing_features": ";".join(sorted(required - features)),
        "missing_signature_features": ";".join(sorted(signature_required - features)),
        "required_refinements": ";".join(sorted(k for k in all_required | signature_required
                                                  if re.fullmatch(r"K\d{6}", k))),
    }
