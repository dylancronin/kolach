# AnnoGuild reference fixtures

Source commit: 291b232907e6568597178e21384ff0f964993a82.

annoguild_parity.json captures 27 deterministic evidence profiles across the
55 refined pathways. Each contains the old evaluator's Boolean calls and
coverage/signature diagnostics. Profiles include empty, complete, focused
markers and random subsets seeded with 1042026.

Initial capture executes the original RuleParser source unchanged, with small
standard-library adapters implementing the exact graph/scalar/lookup operations
it uses. This was necessary because pandas/numpy/networkx could not be installed
in the build environment. The JSON expectations are frozen, not recomputed from
Kolach's new evaluator during tests. To independently verify with real libraries:

```bash
# In an environment containing the original evaluator's dependencies:
python scripts/verify_annoguild_parity.py --annoguild-dir ../AnnoGuild
```

The old get_ko_percent diagnostic raises on an empty option list. Capture records
zero for that diagnostic, consistent with its False Boolean coverage check.
Comammox signature diagnostics expose the old sibling OR-node collision.
New tests compare all Boolean calls and numeric coverage, but deliberately
replace that faulty signature behavior with explicit strict-logic tests.

legacy_calculate_specialisation.py is an unmodified copy of the source classifier.
Tests extract its constants/functions and apply only documented comma, pathway
name and precedence corrections before comparing to Kolach. A dictionary-returning
pd.Series adapter avoids requiring pandas; none of those classifier functions
perform pandas operations. The copy is under the upstream GPL-3.0 license,
included in src/kolach/metabolism/data/LICENSE.GPL-3.0.

The example expected products were generated with the original evaluator and
the corrected original classifier using the supplied marker/CAZy inputs. They
are compared byte-for-byte with Kolach's outputs.
