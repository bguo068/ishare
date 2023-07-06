import pandas as pd
import pybedtools as pb
from subprocess import run
import numpy as np

# %% ============================================
# IBD overlapping analyses via ibdutils compare command


# run ishare ibdutils command
cmd = f"""
cd ..
cargo run -q --release --bin ibdutils  -- compare \
    -g testdata/genome.toml  \
    -s testdata/samples.hap.txt \
    -f tskibd \
    -i testdata/ibd_tskibd \
    -S testdata/samples.hap.txt \
    -F hmmibd \
    -I testdata/ibd_hmmibd \
    -o /tmp/ttttt.csv
"""
r = run(cmd, shell=True, text=True)
if r.returncode != 0:
    print(r.stderr)
    raise Exception("command err")

df1 = pd.read_csv("/tmp/ttttt.csv")
df1


# %%  ===============================================================
# IBD overlapping analyses using pybedtools


lst1, lst2 = [], []

for chrno in range(1, 4):
    # Id1 Id2 Start End Ancestor Tmrca HasMutation
    df1 = pd.read_csv(f"../testdata/ibd_tskibd/{chrno}.ibd", sep="\t")
    df1["Chrom"] = chrno
    lst1.append(df1[["Id1", "Id2", "Chrom", "Start", "End"]])

    # sample1 sample2 chr start end different Nsnp
    df2 = pd.read_csv(f"../testdata/ibd_hmmibd2/{chrno}.hmm.txt", sep="\t")
    df2 = df2[
        (df2["different"] == 0) & (df2["end"] - df2["start"] >= 15000 * 2.0)
    ].iloc[:, :5]
    df2.columns = ["Id1", "Id2", "Chrom", "Start", "End"]
    lst2.append(df2)


# %%
ibd1 = pd.concat(lst1, axis=0)
ibd2 = pd.concat(lst2, axis=0)


# compare results
id1 = ibd1[["Id1", "Id2"]].max(axis=1)
id2 = ibd1[["Id1", "Id2"]].min(axis=1)
bed1 = pd.DataFrame(
    {
        "Chrom": 100000000 * 1 + 100000 * id1 + 100 * id2 + ibd1.Chrom,
        "Start": ibd1["Start"],
        "End": ibd1["End"],
    }
).sort_values(["Chrom", "Start", "End"])
bed1 = pb.BedTool.from_dataframe(bed1)

id1 = ibd2[["Id1", "Id2"]].max(axis=1)
id2 = ibd2[["Id1", "Id2"]].min(axis=1)
bed2 = pd.DataFrame(
    {
        "Chrom": 100000000 * 1 + 100000 * id1 + 100 * id2 + ibd2.Chrom,
        "Start": ibd2["Start"],
        "End": ibd2["End"],
    }
).sort_values(["Chrom", "Start", "End"])
bed2 = pb.BedTool.from_dataframe(bed2)


intersect_1_by_2 = bed1.intersect(bed2, sorted=True, wao=True).to_dataframe()
intersect_1_by_2 = intersect_1_by_2.iloc[:, [0, 1, 2, 6]].copy()
intersect_1_by_2.columns = ["Chrom", "Start", "End", "Bp"]
# summarize
intersect_1_by_2 = (
    intersect_1_by_2.groupby(["Chrom", "Start", "End"])["Bp"].sum().reset_index()
)
intersect_1_by_2["Len"] = intersect_1_by_2["End"] - intersect_1_by_2["Start"]
intersect_1_by_2["Ov"] = intersect_1_by_2["Bp"] / intersect_1_by_2["Len"]
intersect_1_by_2["Cm"] = intersect_1_by_2["Len"] / 15000
intersect_1_by_2 = intersect_1_by_2[["Cm", "Ov"]]
# categorize
intersect_1_by_2["Cm"] = pd.cut(
    intersect_1_by_2["Cm"], [3, 4, 6, 10, 18, np.inf], right=False
)
x = intersect_1_by_2.dropna().groupby("Cm")["Ov"].mean().rename("Ov1By2")

intersect_2_by_1 = bed2.intersect(bed1, sorted=True, wao=True).to_dataframe()
intersect_2_by_1 = intersect_2_by_1.iloc[:, [0, 1, 2, 6]].copy()
intersect_2_by_1.columns = ["Chrom", "Start", "End", "Bp"]
# summarize
intersect_2_by_1 = (
    intersect_2_by_1.groupby(["Chrom", "Start", "End"])["Bp"].sum().reset_index()
)
intersect_2_by_1["Len"] = intersect_2_by_1["End"] - intersect_2_by_1["Start"]
intersect_2_by_1["Ov"] = intersect_2_by_1["Bp"] / intersect_2_by_1["Len"]
intersect_2_by_1["Cm"] = intersect_2_by_1["Len"] / 15000
intersect_2_by_1 = intersect_2_by_1[["Cm", "Ov"]]
# categorize
intersect_2_by_1["Cm"] = pd.cut(
    intersect_2_by_1["Cm"], [3, 4, 6, 10, 18, np.inf], right=False
)
y = intersect_2_by_1.dropna().groupby("Cm")["Ov"].mean().rename("Ov2By1")

df2 = pd.concat([x, y], axis=1)

# %% compare results of two methods

diff1 = np.abs(df1.iloc[:5, 1].values - df2.iloc[:5, 0].values)
diff2 = np.abs(df1.iloc[:5, 2].values - df2.iloc[:5, 1].values)
assert np.all(diff1 < 0.001)
assert np.all(diff2 < 0.001)


# %%
