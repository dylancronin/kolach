# AnnoGuild-derived metabolic rules and classification

The feature adapts definitions and classification logic from
[dylancronin/AnnoGuild](https://github.com/dylancronin/AnnoGuild), pinned at
`291b232907e6568597178e21384ff0f964993a82`.

Sources: resources/custom_input_modules/EMERGE_pathways_module.tsv,
resources/pathways_refined.tsv, scripts/calculate_specialisation.py,
and scripts/paths_graph_parser.py (reference scoring behavior).

AnnoGuild distributes this material under GNU GPL version 3. Its license is
included in src/kolach/metabolism/data/LICENSE.GPL-3.0. The adapted metabolism
code, rules and copied test reference are provided under GPL-3.0. Existing
Kolach files retain their previous licensing status.

The evaluator is rewritten around expressions rather than graphs. Reaction
data have explicit IDs. Corrections and preserved assumptions are recorded
in CHANGELOG.md and docs/specialisation.md.

# DRAM2 CAZy substrate rules

`rules_common.tsv`, `rules_gut.tsv`, `rules_aquatic.tsv` and `rules_eng_sys.tsv`
in `src/kolach/metabolism/data/dram_rules/` are copied unchanged from
[BortonWrightonLabs/dram-rules](https://github.com/BortonWrightonLabs/dram-rules)
at `3257b8604ab6e40c00e811d6ba0211b44dcf5648` (`src/dram_rules/data/`).
They are distributed under the MIT
licence included as `src/kolach/metabolism/data/dram_rules/LICENSE.MIT`.
Kolach evaluates only its CAZy family/EC substrate rules.

# GraftM and GraftM packages

GraftM (GPL-3.0-or-later) is a runtime dependency installed in the
environment; no GraftM code is copied into Kolach.

The GraftM package bundle downloaded by `kolach download --databases graftm`
contains eight packages built for the EMERGE temporal paper
(`temporal_paper_202505/GraftM-analysis`) and the public GraftM package
`7.60.mmoX_sMMO` (https://data.ace.uq.edu.au/public/graftm/7/). The EMERGE
packages' reference sequences include UniProt UniRef90 (release 2022_01)
sequences, available under CC BY 4.0 (https://www.uniprot.org/help/license);
UniProt Consortium, "UniProt: the Universal Protein Knowledgebase". The placement
lookup in `src/kolach/metabolism/data/graftm_placements.tsv` lists every clade
of those packages and is derived from the temporal paper's
`update_annotations.R`, with the corrections described in
docs/specialisation.md.
