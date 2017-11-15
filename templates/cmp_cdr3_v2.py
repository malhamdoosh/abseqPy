import sys
import pandas as pd


cols = ["cdr1", "cdr2", "fr1", "fr2", "fr3", "fr4"]


def main():
    df = pd.read_hdf(sys.argv[1])

    # obtain all unique cdr3s
    cdr3s = set(df['cdr3'].tolist())
    print("Unique cdr3 counts", len(cdr3s))


    with open("v-cdr3-comp-no-strict.csv", "w") as fp:
        # write csv header
        fp.write("cdr3,count," + ','.join(cols) + "\n")

        # for each unqiue cdr3, find all rows(reads) that have this same CDR-H3
        for cdr3 in cdr3s:
            rows = df[df['cdr3'] == cdr3]
            fp.write(cdr3 + "," + str(len(rows)) + "," + ','.join(analyze(rows)) + '\n')
    return 0


def analyze(rows):
    """
    :param rows: all rows that share the SAME CDR3
    :return: number each of unique [cdr1, cdr2, fr1, fr2, fr3, fr4] region counts
    """
    lst = []
    for region in cols:
        lst.append(str(len(set(rows[region].tolist()))))
    return lst


if __name__ == '__main__':
    sys.exit(main())
