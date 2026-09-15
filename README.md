# kolach

Genome annotation targeting KEGG Orthologies (KOs).

`kolach` provides unified workflows to download databases and assign KO identifiers to protein sequences using three complementary methods:
- **KOfam**: Profile HMM searches with bitscore thresholds and heuristic rescue.
- **DeepKOALA**: Deep learning GRU models for complete or fragmented sequences.
- **eggNOG**: Orthology assignment via eggNOG orthologous groups and DIAMOND/MMseqs2.

---

## Installation

```bash
conda env create -f environment.yaml
conda activate kolach
```

---

## Usage

### 1. Download Databases

```bash
kolach download --database-dir /path/to/databases --databases kofam deepkoala eggnog
```

Optional download flags:
- `--deepkoala-release`: Model release tag (`latest` or release tag like `202608`, default: `latest`).

### 2. Annotate Proteins

```bash
kolach annotate \
    --protein-fasta proteins.faa \
    --database-dir /path/to/databases \
    --databases kofam deepkoala eggnog \
    --output-dir kolach_output \
    --threads 8
```

Key options:
- `--deepkoala-model`: `full` (complete proteins, default) or `frag` (metagenomic fragments).
- `--device`: Compute device for DeepKOALA (`auto`, `cpu`, or `cuda`).
- `--batch-size`: Batch size for DeepKOALA inference (default: `64`).
- `--eggnog-mode`: Search mode (`diamond` [default] or `mmseqs`).
- `--eggnog-sensmode`: DIAMOND sensitivity mode (default: upstream sensitive iterative search).
- `--eggnog-min-bitscore` / `--eggnog-max-evalue`: eggNOG retention thresholds (default: `60.0` / `1e-5`).
- `--eggnog-filter-multi`: Resolve multi-KO orthology groups (`disambiguate` [default], `strict`, or `none`).
- `--conflict-strategy`: Adjudication for disjoint tool calls (`multiple` [default] or `priority`).

### 3. Standalone Table Integration

Pre-computed annotation tables can be integrated directly:

```bash
kolach integrate \
    --protein-fasta proteins.faa \
    --kofam-table kolach_output/kofam_annotations.tsv \
    --deepkoala-table kolach_output/deepkoala_annotations.tsv \
    --eggnog-table kolach_output/eggnog_annotations.tsv \
    --output-file kolach_annotations.tsv \
    --evidence-file kolach_evidence.tsv \
    --database-dir /path/to/databases
```

---

## Consensus & Integration Workflow

`kolach` reconciles multi-tool outputs into `{output_dir}/kolach_annotations.tsv` (gene-level summary) and `{output_dir}/kolach_evidence.tsv` (long-form provenance).

### Integration Logic

1. **Confidence Gating & Disambiguation**: Hits must meet method thresholds to count as confident votes (KOfam bit score $\ge$ profile threshold, DeepKOALA probability $\ge$ threshold, eggNOG bit score $\ge$ 60 and E-value $\le$ 1e-5). Sub-threshold predictions are tracked as candidates. eggNOG multi-KO hits are disambiguated against other tools' predictions.
2. **Consensus by Confident Call Count (3, 2, 1, 0)**:
   - **3 Confident Calls**: All 3 agree $\rightarrow$ `unanimous`. Exactly 2 agree $\rightarrow$ `majority` (the 3rd unselected call is preserved in `alternative_kos`). All 3 differ $\rightarrow$ disjoint conflict.
   - **2 Confident Calls**: Both agree $\rightarrow$ `majority` (or `unanimous` if only 2 tools were run). Disagree $\rightarrow$ disjoint conflict.
   - **1 Confident Call**: Supported by another tool's sub-threshold candidate $\rightarrow$ `single_tool_with_candidate`; sole call $\rightarrow$ `single_tool`.
   - **0 Confident Calls**: 2 independent sub-threshold candidates agree $\rightarrow$ `orthogonal_dual_candidate`. KOfam heuristic rescue $\rightarrow$ `single_tool` (`evidence = kofam(rescued)`). Otherwise $\rightarrow$ `unannotated`.
3. **Conflict Resolution (`--conflict-strategy`)**:
   - `multiple` *(default)*: Assigns `consensus_level = conflict`, sets `accepted_ko = '-'`, and records all conflicting KOs in `alternative_kos`.
   - `priority`: Assigns `consensus_level = conflict_priority`, selecting the top tool's KO (`kofam > eggnog > deepkoala`) for `accepted_ko` and moving unselected KOs to `alternative_kos`.
4. **Strict Single-KO Policy**: Downstream pathway tools require single-KO calls. `accepted_ko` strictly contains exactly one KO (or `-` if unannotated, multi-KO, or conflicting). All minority, dropped, or conflicting KOs are retained in `alternative_kos`.

### Workflow Diagram

```mermaid
flowchart TD
    Start["Gene Predictions<br/>(KOfam, DeepKOALA, eggNOG)"] --> Count{"Confident calls<br/>meeting thresholds?"}

    %% 3 calls
    Count -->|"3 calls"| C3{"Agreement?"}
    C3 -->|"3 agree"| Unanimous["consensus_level: unanimous"]
    C3 -->|"2 agree"| Majority["consensus_level: majority<br/>(3rd KO → alternative_kos)"]
    C3 -->|"0 agree"| Conflict{"--conflict-strategy"}

    %% 2 calls
    Count -->|"2 calls"| C2{"Agreement?"}
    C2 -->|"2 agree"| Maj2["consensus_level: majority<br/>(unanimous if 2 tools run)"]
    C2 -->|"0 agree"| Conflict

    %% Conflict resolution
    Conflict -->|"multiple (default)"| ConfM["consensus_level: conflict<br/>(accepted_ko = '-', all to alternative_kos)"]
    Conflict -->|"priority"| ConfP["consensus_level: conflict_priority<br/>(accept top tool: kofam > eggnog > deepkoala)"]

    %% 1 call
    Count -->|"1 call"| C1{"Candidate support<br/>from other tool?"}
    C1 -->|"Yes"| SingleCand["consensus_level: single_tool_with_candidate"]
    C1 -->|"No"| Single["consensus_level: single_tool"]

    %% 0 calls
    Count -->|"0 calls"| C0{"Sub-threshold / Rescue?"}
    C0 -->|"2 candidates agree"| DualCand["consensus_level: orthogonal_dual_candidate"]
    C0 -->|"KOfam rescued"| Rescued["consensus_level: single_tool<br/>(evidence: kofam(rescued))"]
    C0 -->|"None"| Unannotated["consensus_level: unannotated<br/>(accepted_ko = '-')"]
```

### Consensus Categories (`consensus_level`)

| Level | Description | Output Routing |
| :--- | :--- | :--- |
| `unanimous` | All active methods confidently agree on the KO. | `accepted_ko` = single KO |
| `majority` | At least 2 methods confidently agree on the KO. | `accepted_ko` = agreed KO (`alternative_kos` = unselected 3rd KO) |
| `single_tool_with_candidate` | 1 confident call corroborated by another tool's candidate. | `accepted_ko` = single KO |
| `single_tool` | Exactly 1 confident method (or KOfam rescue) without corroboration. | `accepted_ko` = single KO |
| `orthogonal_dual_candidate` | Sub-threshold candidates from 2 methods agree on the same KO. | `accepted_ko` = agreed KO |
| `conflict` | Disjoint calls under `multiple` (default). | `accepted_ko` = `-` (`alternative_kos` = all conflicting KOs) |
| `conflict_priority` | Disjoint calls resolved by method hierarchy. | `accepted_ko` = top tool KO (`alternative_kos` = unselected KOs) |
| `unannotated` | No method produced a confident or corroborated assignment. | `accepted_ko` = `-` (`alternative_kos` = `-`) |

---

## Output Schemas

### 1. Gene Summary Table (`kolach_annotations.tsv`)

| Column | Description |
| :--- | :--- |
| `gene_id` | Protein identifier from FASTA (preserves input order). |
| `accepted_ko` | Accepted single KEGG Orthology identifier (`-` when unannotated, conflicting, or multi-KO). |
| `alternative_kos` | Unresolved multi-KO or conflicting candidate KOs. |
| `definition` | Official KEGG functional definition for `accepted_ko`. |
| `alternative_definition` | Functional definitions corresponding to `alternative_kos`. |
| `consensus_level` | Classification status (`unanimous`, `majority`, `single_tool`, `conflict`, etc.). |
| `evidence` | Methods and candidates supporting the assignment. |
| `kofam_ko` | Confident KO(s) assigned by KOfam (or `-`). |
| `kofam_score_type` | Applicable profile score type (`full` or `domain`). |
| `kofam_bit_score` | KOfam bit score (selected by score type). |
| `kofam_evalue` | KOfam E-value (selected by score type). |
| `kofam_assignment` | KOfam assignment type (`threshold` or `rescued`). |
| `kofam_threshold` | KOfam bit score threshold. |
| `deepkoala_ko` | Confident KO meeting DeepKOALA threshold (or `-`). |
| `deepkoala_candidate_ko` | Raw candidate KO predicted by DeepKOALA. |
| `deepkoala_score` | DeepKOALA model prediction probability. |
| `deepkoala_threshold` | DeepKOALA model confidence threshold. |
| `eggnog_ko` | Confident filtered KO(s) meeting score and E-value cutoffs (or `-`). |
| `eggnog_candidate_ko` | Raw candidate KO(s) predicted by eggNOG. |
| `eggnog_bit_score` | eggNOG alignment bit score. |
| `eggnog_evalue` | eggNOG alignment E-value. |

### 2. Evidence Table (`kolach_evidence.tsv`)

Long-form table with one row per gene $\times$ KO $\times$ method, preserving full provenance:

| Column | Description |
| :--- | :--- |
| `gene_id` | Protein identifier. |
| `ko` | Specific KEGG Orthology identifier. |
| `method` | Annotation method (`kofam`, `deepkoala`, or `eggnog`). |
| `bit_score` | Full-sequence alignment bit score. |
| `e_value` | Full-sequence alignment E-value. |
| `domain_bit_score` | Domain bit score for HMM methods. |
| `domain_e_value` | Domain E-value for HMM methods. |
| `score_type` | Threshold score type (`full`, `domain`, or `seed_hit`). |
| `bit_score_threshold` | Bit score cutoff for method (`0.75 * threshold` for rescued KOfam, profile cutoff for primary KOfam, min bitscore for eggNOG, `-` for DeepKOALA). |
| `e_value_threshold` | E-value cutoff for method (`1e-5` for rescued KOfam and eggNOG, `-` for primary KOfam and DeepKOALA). |
| `deepkoala_probability` | DeepKOALA prediction probability score. |
| `deepkoala_threshold` | DeepKOALA confidence threshold. |
| `eggnog_shared_seed_hit` | Boolean indicating whether eggNOG alignment metrics are shared across multiple candidate KOs (`True` for multi-KO eggNOG rows, `False` otherwise). |
| `original_status` | Pre-consensus status (`threshold_passing`, `heuristic_rescued`, `below_threshold`, `unknown`). |
| `cross_method_status` | Integration outcome (`accepted`, `heuristic_rescued`, `rescued`, `disambiguated_retained`, `disambiguated_dropped`, `alternative`, `conflict`, `unresolved_multi_ko`, `below_threshold_unrescued`, `unselected`). |
| `consensus_level` | Overall consensus classification for the gene (`unanimous`, `majority`, `single_tool`, etc.). |
| `supporting_methods` | Other methods corroborating this specific KO. |
