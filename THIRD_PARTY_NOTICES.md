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
# CAZy reference definitions

`src/kolach/metabolism/data/cazy_definitions.tsv` is an unchanged copy of
`resources/DRAM_CAZy_definitions.tsv` from dylancronin/AnnoGuild commit
`291b232907e6568597178e21384ff0f964993a82` (Git blob
`508a9bb2e6502697739bef0e6225e36144953d49`). The existing bundled GPL-3.0
license applies to these AnnoGuild assets. dbCAN is an optional external
dependency, pinned to its released version 5.2.9; its annotation algorithms
are not copied into Kolach.
