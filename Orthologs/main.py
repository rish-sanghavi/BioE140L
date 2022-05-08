from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
Entrez.email = "" #FIXME, add in your own email.
import os
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

def setUp():
    organisms = []
    taxida = []
    if os.path.isdir("orthologs") == False:
        os.makedirs("orthologs")
    while (True):
        k = input("Please enter in all the names of the organisms you would like to find orthologs for. Type done when finished.")
        if k == "done" or k == "Done":
            break
        else:
            taxidtemp = input("Please enter the taxid of the organism " + k + ".")
            taxida.append(taxidtemp)
            organisms.append(k)
    findOrthologs(organisms, taxida)


def findOrthologs(organisms, taxida):
    seq_input = input("Please enter the protein sequence you would like to find Orthologs for.")
    print("Finding hits! This will take a few minutes so please be patient =)")
    count = 0
    for x in organisms:
        try:
            taxid = taxida[count]
            result = NCBIWWW.qblast("blastp", "nr", seq_input, entrez_query=("txid" + taxid + "[ORGN]"))
            savefile = open(x + "cdsp.xml", "w")
            savefile.write(result.read())
            savefile.close()
            blastres = open(x + "cdsp.xml", "r")
            blastrecord = NCBIXML.read(blastres)
            proteinid = blastrecord.alignments[0].hit_id
            count2 = 0
            realid2 = ""
            for u in proteinid:
                if (u == "|"):
                    count2 += 1
                elif (count2 == 1):
                    realid2 = realid2 + u
            proteinseq = Entrez.efetch(db="protein", id=realid2, retmode="xml")
            hitseq = Entrez.read(proteinseq)[0]["GBSeq_sequence"]
            result2 = NCBIWWW.qblast("tblastn", "nt", hitseq, entrez_query=("txid" + taxid + "[ORGN]"))
            savefile2 = open(x + "cdstn.xml", "w")
            savefile2.write(result2.read())
            savefile2.close()
            blastres2 = open(x + "cdstn.xml", "r")
            blastrecord2 = NCBIXML.read(blastres2)
            accession = blastrecord2.alignments[0].accession
            id = blastrecord2.alignments[0].hit_id
            start = blastrecord2.alignments[0].hsps[0].sbjct_start
            end = blastrecord2.alignments[0].hsps[0].sbjct_end
            realid = ""
            count3 = 0
            for f in id:
                if (f == "|"):
                    count3 += 1
                elif (count3 == 1):
                    realid = realid + f
            genomelength = blastrecord2.alignments[0].length
            if (start < 200):
                start2 = 0
            else:
                start2 = start - 200
            if (end + 200 > genomelength):
                end2 = genomelength
            else:
                end2 = end + 200
            nukeseq = Entrez.efetch(db="nucleotide", id=realid, retmode="xml", seq_start=start2, seq_stop=end2)
            genomeseq = Entrez.read(nukeseq)[0]["GBSeq_sequence"]
            data = dataProcess(genomeseq, start, start2, end, accession)
            writeData(data, x)
            count += 1
        except:
            print("Could not find Ortholog for " + x + ".")
            count += 1
    print("Finsihed finding Orthologs!")


def dataProcess(data, start, start2, end, accession):
    startcodons = ["atg", "ttg", "ctg", "gtg"]
    stopcodons = ["tag", "tga", "taa"]
    preseqstart = 0
    preseqend = (start - start2)
    cdsstart = preseqend
    cdsend = (end - start) + cdsstart + 1
    endseqstart = cdsend
    endseqend = len(data) + 1
    pre = Seq((data[preseqstart:preseqend]))
    cds = Seq((data[cdsstart:cdsend]))
    post = Seq((data[endseqstart:endseqend]))
    if (cds[0:3] in startcodons):
        if (cds[len(cds) - 3:len(cds)] in stopcodons):
            return [str(pre), str(post), str(cds), accession]
        else:
            residue = ""
            codon = ""
            found = False
            for i in post:
                if (len(codon) == 3):
                    if (codon in stopcodons):
                        found = True
                        break
                    else:
                        codon = ""
                residue = residue + i
                codon = codon + i
            if not found:
                return
            else:
                cds = cds + residue
                post = post[len(residue):]
                return [str(pre), str(post), str(cds), accession]
    elif (cds.reverse_complement()[0:3] in startcodons):
        cds = cds.reverse_complement()
        post2 = post
        post = pre.reverse_complement()
        pre = post2.reverse_complement()
        if (cds[len(cds) - 3:len(cds)] in stopcodons):
            return [str(pre), str(post), str(cds), accession]
        else:
            residue = ""
            codon = ""
            found = False
            for i in post:
                if (len(codon) == 3):
                    if (codon in stopcodons):
                        found = True
                        break
                    else:
                        codon = ""
                residue = residue + i
                codon = codon + i
            if not found:
                return
            else:
                cds = cds + residue
                post = post[len(residue):]
                return [str(pre), str(post), str(cds), accession]
    else:
        return

def writeData(data, x):
    preseq = data[0]
    postseq = data[1]
    cds = data[2]
    accession = data[3]
    open("orthologs/" + x + ".txt", "w")
    with open("orthologs/" + x + ".txt", "w") as f:
        f.write("Accession#, Presequence, CDS, Postsequence" + "\n")
        f.write(accession + "\t")
        f.write(preseq + "\t")
        f.write(cds + "\t")
        f.write(postseq + "\t")
        f.write("\n")
    print("Finished finding ortholog for " + x + ".")

setUp()

