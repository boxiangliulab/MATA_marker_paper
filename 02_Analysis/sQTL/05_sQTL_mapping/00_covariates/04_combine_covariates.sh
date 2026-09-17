#!/bin/bash
# Assemble the covariate files.
#
# Combine the three covariate sources into the single file that ISSAC reads
# for each cell type and depot. This happens in three passes:
#
#   (a) per donor : genotype PCs (step 01) + age/sex/diabetes/BMI (step 02)
#   (b) per sample: the result of (a) attached to every meta-cell of that
#                   donor, next to the splicing PCs (step 03). Meta-cell
#                   names are <donor>:<index>, so the donor is the part
#                   before the colon.
#   (c) formatting: the table is transposed (ISSAC expects one column per
#                   sample) and the categorical covariates are encoded
#                   numerically: sex M/F -> 1/0, diabetes T2DM/x -> 1/0.

set -euo pipefail
cd "$(dirname "$0")"

# Paths
MODEL_DIR=../../data/05_sQTL_mapping

pc_dir="${MODEL_DIR}/PCs"
p=5   # splicing PCs

# --- (a) Donor level: genotype PCs + the other donor covariates --------------
for g in 10 6; do
    gPCs="${pc_dir}/genotype_PCs.${g}.txt"
    merge_PCs="${pc_dir}/merge_PCs.${g}.txt"

    awk -v G="${g}" '
    BEGIN { FS = OFS = "\t" }
    FNR == 1 && NR == FNR { next }          # skip the genotype-PC header
    NR == FNR { a[$1] = $0; next }          # genotype PCs, keyed by donor
    FNR == 1 {                              # covariate header + gPC columns
        header = $0
        for (i = 1; i <= G; i++) header = header "\tgPC" i
        print header
        next
    }
    $1 in a {                               # covariates + that donor gPCs
        split(a[$1], b, "\t")
        line = $0
        for (i = 2; i <= G + 1; i++) line = line OFS b[i]
        print line
    }' "${gPCs}" "${pc_dir}/sample_covariates.txt" > "${merge_PCs}"

    echo "Wrote ${merge_PCs}"
done

# --- (b) + (c) Sample level, per cell type and depot -------------------------
for depot in SAT VAT; do
    echo "=== ${depot} ==="
    if [ "${depot}" == "SAT" ]; then g=10; else g=6; fi

    cut -f1 "../../reference/celltypes.QTL.${depot}.txt" | while read -r ct; do
        merge_PCs="${pc_dir}/merge_PCs.${g}.txt"
        splicing_PCs="${pc_dir}/${ct}/${depot}/splicing_PCs.${p}.txt"
        all_PCs="${pc_dir}/${ct}/${depot}/all_PCs.g${g}.p${p}.txt"
        all_PCs_format="${pc_dir}/${ct}/${depot}/all_PCs.g${g}.p${p}.format.txt"

        if [ ! -s "${splicing_PCs}" ]; then
            echo "  skip ${ct}: ${splicing_PCs} not found"
            continue
        fi

        echo "--- ${ct} ---"
        mkdir -p "${pc_dir}/${ct}/${depot}"

        # (b) attach the donor-level covariates to every meta-cell
        awk '
        BEGIN { FS = OFS = "\t" }
        NR == FNR {                          # donor-level table
            if (FNR == 1) { header2 = $0 } else { a[$1] = $0 }
            next
        }
        FNR == 1 {                           # splicing-PC header + donor columns
            split(header2, h2, "\t")
            for (i = 1; i <= NF; i++) printf "%s\t", $i
            for (i = 2; i <= length(h2); i++) {
                printf "%s", h2[i]
                if (i < length(h2)) printf "\t"
            }
            print ""
            next
        }
        {                                    # <donor>:<index> -> donor
            split($1, id, ":")
            if (id[1] in a) {
                split(a[id[1]], b, "\t")
                printf "%s", $0
                for (i = 2; i <= length(b); i++) printf "\t%s", b[i]
                print ""
            }
        }' "${merge_PCs}" "${splicing_PCs}" > "${all_PCs}"

        # (c) transpose (samples become columns) and encode the categorical
        #     covariates; the header line is dropped because the sample names
        #     are the first column and therefore become the first row.
        tail -n +2 "${all_PCs}" | awk '
        {
            for (i = 1; i <= NF; i++) a[NR, i] = $i
        }
        NF > p { p = NF }
        END {
            for (i = 1; i <= p; i++) {
                str = a[1, i]
                for (j = 2; j <= NR; j++) str = str "\t" a[j, i]
                print str
            }
        }' | sed 's/T2DM/1/g' | sed 's/x/0/g' | sed 's/M/1/g' | sed 's/F/0/g' \
           > "${all_PCs_format}"
    done
done

echo "Covariate files written under ${pc_dir}"
