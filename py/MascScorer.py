import pandas as pd
import string, StringIO, re

class MascScorer:
    """A class to produce standardized t scores from raw masc scores, when combined with gender and age information"""
    
    _tScoreColumn=6

    _femaleMascTable=pd.read_csv(StringIO.StringIO(
        """C8To11YearOlds_MASC,C12To15YearOlds_MASC,C16To19YearOlds_MASC,C8To11YearOlds_ADI,C12To15YearOlds_ADI,C16To19YearOlds_ADI,T
        113-114,102,104,,,,90
        111-112,100-101,102-103,,29,29,89
        110,99,101,,,,88
        108-109,97-98,99-100,,28,28,87
        107,96,98,30,,,86
        105-106,95,96-97,,,,85
        103-104,93-94,95,29,27,27,84
        102,92,94,,,,83
        100-101,90-91,92-93,28,26,26,82
        99,89,91,,,,81
        97-98,87-88,89-90,27,,25,80
        95-96,86,88,,25,,79
        94,84-85,86-87,26,,,78
        92-93,83,85,,24,24,77
        91,82,83-84,,,,76
        89-90,80-81,82,25,23,23,75
        87-88,79,80-81,,,,74
        86,77-78,79,24,,,73
        84-85,76,77-78,,22,22,72
        83,74-75,76,23,,,71
        81-82,73,75,,21,21,70
        79-80,71-72,73-74,22,,,69
        78,70,72,,20,,68
        76-77,69,70-71,21,,20,67
        74-75,67-68,69,,,,66
        73,66,67-68,,19,19,65
        71-72,64-65,66,20,,,64
        70,63,64-65,,18,,63
        68-69,61-62,63,19,,18,62
        66-67,60,61-62,,,,61
        66-67,58-59,60,18,17,17,60
        63-64,57,59,,,,59
        62,56,57-58,17,16,,58
        60-61,54-55,56,,,16,57
        58-59,53,54-55,16,15,,56
        57,51-52,53,,,15,55
        55-56,50,51-52,,,,54
        54,48-49,50,15,14,,53
        52-53,47,48-49,,,14,52
        50-51,45-46,47,14,13,,51
        49,44,45-46,,,13,50
        47-48,43,44,13,,,49
        46,41-42,42-43,,12,,48
        44-35,40,41,12,,12,47
        42-43,38-39,40,,11,,46
        41,37,38-39,11,,11,45
        39-40,35-36,37,,10,,44
        38,34,35-36,10,,,43
        36-37,32-33,34,,,10,42
        34-35,31,32-33,,9,,41
        33,30,31,9,,9,40
        31-32,28-29,29-30,,8,,39
        30,27,28,8,,,38
        28-29,25-26,26-27,,7,8,37
        26-27,24,25,7,,,36
        25,22-23,24,,,7,35
        23-24,21,22-23,6,6,,34
        22,19-20,21,,,,33
        20-21,18,19-20,5,5,6,32
        18-19,17,18,,,,31
        17,15-16,16-17,,,5,30
        15-16,14,15,4,4,,29
        13-14,12-13,13-14,,,,28
        12,11,12,3,3,4,27
        10-11,9-10,10-11,,,,26
        9,8,9,2,2,3,25"""))
    
    _maleMascTable=pd.read_csv(StringIO.StringIO(
        """C8To11YearOlds_MASC,C12To15YearOlds_MASC,C16To19YearOlds_MASC,C8To11YearOlds_ADI,C12To15YearOlds_ADI,C16To19YearOlds_ADI,T
        105-106,94,98,30,,28,90
        104,92-93,96-97,,27,,89
        102-103,91,95,29,,,88
        100-101,90,93-94,,26,27,87
        99,88-99,92,28,,,86
        97-98,87,90-91,,25,26,85
        96,85-86,89,27,,,84
        94-95,84,87-88,,,25,83
        92-93,82-83,86,26,24,,82
        91,81,84-85,,,,81
        89-90,80,83,,23,24,80
        88,78-79,81-82,25,,,79
        86-87,77,80,,,23,78
        85,75-76,78-79,24,22,,77
        83-84,74,77,,,22,76
        81-82,72-73,75-76,23,21,,75
        80,71,74,,,,74
        78-79,69-70,72-73,22,,21,73
        77,68,71,,20,,72
        75-76,67,69-70,21,,20,71
        73-74,65-66,68,,19,,70
        72,64,66-67,,,,69
        70-71,62-63,65,20,,19,68
        67-68,59-60,62,19,,18,66
        66,58,60-61,,17,,65
        64-65,57,59,18,,17,64
        62-63,55-56,57-58,,,,63
        61,54,56,17,16,,62
        59-60,52-53,54-55,,,16,61
        58,51,53,16,15,,60
        56-57,49-50,51-52,,,15,59
        54-55,48,50,,14,,58
        53,47,48-49,15,,,57
        51-52,45-46,47,,,14,56
        50,44,45-46,14,13,,55
        48-49,42-43,44,,,13,54
        47,41,42-43,13,12,,53
        45-46,39-40,41,,,12,52
        43-44,38,39-40,12,,,51
        42,36-37,38,,11,,50
        40-41,35,36-37,11,,11,49
        39,34,35,,10,,48
        37-38,32-33,33-34,,,10,47
        35-36,31,32,10,,,46
        34,29-30,30-31,,9,9,45
        32-33,28,29,9,,,44
        31,26-27,27-28,,8,,43
        29-30,25,26,8,,8,42
        28,24,24-25,,,,41
        26-27,22-23,23,7,7,7,40
        24-25,21,21-22,,,,39
        23,19-20,20,6,6,,38
        21-22,18,18-19,,,6,37
        20,16-17,17,,5,,36
        18-19,15,15-16,5,,5,35
        16-17,14,14,,,,34
        15,12-13,12-13,4,4,4,33
        13-14,11,11,,,,32
        12,9-10,9-10,3,3,,31
        10-11,8,7-8,,,3,30
        9,6-7,6,2,,,29
        7-8,5,4-5,,2,2,28
        5-6,3-4,3,1,,,27
        4,2,1-2,,1,1,26
        <4,<2,0,<1,<1,<1,25"""))
    
    def __init__(self):
        #print self._femaleMascTable.columns
        #print self._maleMascTable.columns        
        pass

    def scoreMasc(self, inGender, inAge, inRawMascScore, inDebug=False):

        if inDebug:
            print "inGender=", inGender, "inAge=", inAge, "inRawMascScore=", inRawMascScore
    
        if inGender is None or inAge is None or inRawMascScore is None:
            return None

        if inAge < 8 or inAge > 19:
            print "*** scoreMasc: Cannot score MASC for a person of age %f: outside bounds (8, 19). Returning None" % inAge
            return None

        if inRawMascScore < 0:
            return None
        
        if string.lower(inGender)[0]=="m":
            if inRawMascScore < 0 or inRawMascScore > 106:
                print "*** scoreMasc: Cannot convert a MASC score of %f to a tscore: outside bounds (0, 106). Returning None." % inRawMascScore
                return None
            return self._scoreMascHelper(self._maleMascTable, inGender, inAge, inRawMascScore)
        elif string.lower(inGender)[0]=="f":
            if inRawMascScore < 8 or inRawMascScore > 114:
                print "*** scoreMasc: Cannot convert a MASC score of %f to a tscore: outside bounds (8, 114). Returning None." % inRawMascScore
                return None
            return self._scoreMascHelper(self._femaleMascTable, inGender, inAge, inRawMascScore)
        else:
            print "*** scoreMasc: Unknown gender: %s. Returning None" % inGender
            return None

    def _scoreMascHelper(self, inTable, inGender, inAge, inRawMascScore):

        tScore=None
    
        if (inAge >= 8 and inAge <= 11):
            tableColumn=0
        elif inAge >= 12 and inAge <= 15:
            tableColumn=1
        elif inAge >= 16 and inAge <= 19:
            tableColumn=2
        else:
            msg="Cannot score MASC for a person of age: %f" % inAge
            raise Error(msg)
    
    ##print(sprintf("Using tableColumn %d\n", tableColumn))

    ## this is extremely slow but we engage in a linear search for the
    ## appropriate row of the table with the T score
    
        for i in xrange(0, inTable.shape[0]):
            ##print (sprintf("i=%d", i))
            cell=inTable.iloc[i, tableColumn]
            ##print(cell)

            if re.search("-", cell):
                ## this is a cell with a range of values in it
                range=str.split(cell, "-")
                if inRawMascScore >= int(range[0]) and inRawMascScore <= int(range[1]):
                    tScore=inTable.iloc[i, self._tScoreColumn]
                    break
            elif re.search("<", cell):
                ## the only relational operation in any cell is < so don't worry
                ## about other relational operators
                end=str.split(cell, "<")[1]
                if inRawMascScore < int(end):
                    tScore=inTable.iloc[i, self._tScoreColumn]
                    break
            else:
                ## it must be a single number on its own
                if inRawMascScore == int(cell):
                    tScore=inTable[i, self._tScoreColumn]
                    break
         ## end of for i in xrange(1, inTable.shape[0]):
        return int(tScore)


if __name__ == "__main__":
    xx=MascScorer()
    print xx.scoreMasc("f", 15, 33)

