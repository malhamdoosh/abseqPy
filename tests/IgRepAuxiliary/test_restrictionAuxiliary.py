from abseqPy.IgRepAuxiliary.restrictionAuxiliary import *


def test_loadRestrictionSites(tmpdir):
    tmp_file = str(tmpdir.mkdir("data").join("restruction_sites.txt"))
    with open(tmp_file, 'w') as fp:
        fp.write("ENZ1\tAAAAA\n"
                 "#comment\n"
                 "ENZ2\tAAAAB\n"
                 "ENZ3\tRYS\n"
                 "ENZ4\tNNN\n"
                 "ENZ5\tXFI\n")
        fp.flush()

    # todo: should this block be inside the context manager, assuming that pytest would destroy the datafile?
    sites = loadRestrictionSites(tmp_file)

    # make sure that the IUPAC letters are translated to ACGT regex form:
    assert "ENZ1" in sites and sites["ENZ1"] == "AAAAA"
    assert "ENZ2" in sites and sites["ENZ2"] == "AAAA[CGT]"
    assert "ENZ3" in sites and sites["ENZ3"] == "[AG][CT][GC]"
    assert "ENZ4" in sites and sites["ENZ4"] == "..."

    # make sure that unknown IUPAC letters are left as-is
    assert "ENZ5" in sites and sites["ENZ5"] == "XFI"


def test_replaceIUPACLetters():
    # normal usages
    assert replaceIUPACLetters("AAAA") == "AAAA"
    assert replaceIUPACLetters("NNVH") == "..[ACG][ACT]"
    assert replaceIUPACLetters("ACGT") == "ACGT"

    # when letters are not found in the IUPAC mapping, just use it as it is
    assert replaceIUPACLetters("01239") == "01239"

    # works when letters are lowercase too
    assert replaceIUPACLetters("AAAA") == replaceIUPACLetters("aaaa")
    assert replaceIUPACLetters("NNVH") == replaceIUPACLetters("nnvh")
    assert replaceIUPACLetters("ACGT") == replaceIUPACLetters("acgt")
