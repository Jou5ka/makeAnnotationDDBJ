# -*- coding: utf-8 -*-

"""
USAGE
python annotation.py gff function fasta rrnammer > annotation

"""

import sys
from Bio import SeqIO
from itertools import groupby
import re

def header():

	print("""COMMON\tDATATYPE\t\ttype\tWGS
\tKEYWORD\t\tkeyword\tWGS
\t\t\tkeyword\tSTANDARD_DRAFT
\tDBLINK\t\tproject\tPRJDB8156
\t\t\tbiosample\tSAMD00166938
\tSUBMITTER\t\tab_name\tTakahashi,H.
\t\t\tab_name\tYaguchi,T.
\t\t\tcontact\tHiroki Takahashi
\t\t\temail\thiroki.takahashi@chiba-u.jp
\t\t\tphone\t81-43-226-2822
\t\t\tfax\t81-43-226-2486
\t\t\tinstitute\tChiba University
\t\t\tdepartment\tMedical Mycology Research Center
\t\t\tcountry\tJapan
\t\t\tstate\tChiba
\t\t\tcity\tChiba
\t\t\tstreet\t1-8-1 Inohana
\t\t\tzip\t260-8673
\tREFERENCE\t\ttitle\tDraft genome sequence of Aspergillus udagawae IFM 46972
\t\t\tab_name\tTakahashi,H.
\t\t\tab_name\tYaguchi,T.
\t\t\tstatus\tUnpublished
\t\t\tyear\t2019
\tST_COMMENT\t\ttagset_id\tGenome-Assembly-Data
\t\t\tAssembly Method\tSPAdes v. 3.12.0
\t\t\tGenome Coverage\t38x
\t\t\tSequencing Technology\tHiSeq 1500
\tDATE\t\thold_date\t20190630""")

def codonCheck(startEndList, seq, strand):
	STOPCODON=["TAA","TAG","TGA"]

	num=len(startEndList)-1
	transcript=""

	if strand=="+":
		for item in startEndList:
			start=int(item[0])-1
			end=int(item[1])
			transcript=transcript+seq[start:end]

		sLen=len(transcript)
		startCodon=transcript[0:3]
		stopCodon=transcript[(sLen-3):sLen]

		if startCodon!="ATG":
			startEndList[0][0]="<"+startEndList[0][0]

		if stopCodon not in STOPCODON:
			startEndList[num][1]=">"+startEndList[num][1]


	if strand=="-":
		for item in startEndList:
			start=int(item[0])-1
			end=int(item[1])
			tmp=seq[start:end]
			transcript=transcript+tmp.reverse_complement()

		sLen=len(transcript)
		startCodon=transcript[0:3]
		stopCodon=transcript[(sLen-3):sLen]

		if startCodon!="ATG":
			startEndList[0][1]=">"+startEndList[0][1]

		if stopCodon not in STOPCODON:
			startEndList[num][0]="<"+startEndList[num][0]

	return startEndList

def gettRNA(file):	#file: gff

	tRNA={}
	## get tRNA ID
	with open(file) as f:
		for line in f:
			line=line.rstrip()
			s=line.split()

			try:
				if s[2]=="tRNA":
					feat=s[8].split(";")[0]
					id=feat.split("ID=")[1]
					tRNA[id]=s[0]
			except IndexError:
				pass
	#tRNA: {IFM46972_06391-RA:NODE_43_length_239763_cov_15.224840

	## get tRNA product
	with open(file) as f:
		for line in f:
			line=line.rstrip()
			s=line.split()

			try:
				if s[2]=="gene":
					feat=s[8].split(";")[0]
					id=feat.split("ID=")[1]

					key=id+"-RA"		#IFM46972_06391-RA
					if key in tRNA:
						tmp=s[8].split(";")[2] #Alias=trnascan-NODE_43_length_239763_cov_15.224840-noncoding-Gly_GCC-gene-0.92
						tmp2=tmp.split("-noncoding-")[1] #Gly_GCC-gene-0.92
						product="tRNA-"+tmp2.split("_")[0] #tRNA-Gly
						tRNA[key]=product
			except IndexError:
				pass

	##
	tRNADict={}
	with open(file) as f:
		for line in f:
			line=line.rstrip()
			s=line.split()

			try:
				if s[2]=="exon":
					feat=s[8].split(":")[0]
					id=feat.split("ID=")[1]
					if id in tRNA:
						product=tRNA[id]
						key=id+"|"+s[6]+"|"+s[0]+"|"+product
						val=s[3]+".."+s[4]
						if key in tRNADict:
							tRNADict[key]+="|"+val
						else:
							tRNADict[key]=val

			except IndexError:
				pass

	##
	tRNA_data={}
	for k, v in tRNADict.items():
		id,strand,scf,product=k.split("|")
		s=v.split("|")
		pos=""

		if strand=="+":
			if len(s)==1:
				pos=s[0]
			else:
				pos="join("+",".join(s)+")"

		else:
			s.reverse()
			if len(s)==1:
				pos=s[0]
			else:
				pos="complement(join("+",".join(s)+"))"

		val="\ttRNA\t"+pos+"\tproduct\t"+product+"\n\t\t\tlocus_tag\t"+id.split("-")[0]

		if scf in tRNA_data:
			tRNA_data[scf]+="|"+val
		else:
			tRNA_data[scf]=val
	return tRNA_data

def getProt(file):
	## ID	gi	annotation	species
	## IFM46972_11423-RA       gi|849276611|dbj|GAO81050.1|    uncharacterized aminotransferase C1771.03c      Aspergillus udagawae
	## IFM46972_11426-RA       gi|849276608|dbj|GAO81047.1|    galactose oxidase       Aspergillus udagawae
	## IFM46972_11424-RA       gi|849276610|dbj|GAO81049.1|    putative branched-chain-amino-acid aminotransferase     Aspergillus udagawae

	unknowList=("SD08430p", "Similar to Eukaryotic translation initiation factor 5A; acc. no. O94083")
	lowercase=("Acylphosphatase","Ankyrin","Choline","Fibronectin","Fungal","Meiotic","Metallo","Mitochondria","Neopullulanase","Patatin","Phospho","Protein","Putative","Transcription","Trihydro","Uncharacterized")

	f=open(file, 'Ur')
	a={}
	for line in f:
		line=line.rstrip()
		b=line.split("\t")
		if not b[2].startswith("hypothetical protein"):
			if b[2] in unknowList:
				b[2]="hypothetical protein"
			if b[2].startswith(lowercase):
				b[2]=b[2][0].lower()+b[2][1:]
			b[2]=b[2].replace("8s rRNA","8S rRNA")
			a[b[0]]=b[2]
		else:
			a[b[0]]="hypothetical protein"
	f.close()
	return a

def getGap(file):
	gapData={}
	with open(file) as f:
		for record in SeqIO.parse(f, "fasta"):
			for match in re.finditer('N+', str(record.seq)):
				out=str(match.start()+1)+".."+str(match.end())
				if record.id in gapData:
					gapData[record.id]+="|"+out
				else:
					gapData[record.id]=out

	return gapData

def getrRNA(file):

	rRNAdata={}

	with open(file) as f:
		for line in f:
			line=line.rstrip()

			if line.startswith("#"):
				continue

			s=line.split("\t")
				
			## replace "S"
			s[8]=s[8].replace("s", "S")
			prod=s[8].split("_")[0]+" rRNA"
			pos=""
			if s[6]=="+":
				pos=s[3]+".."+s[4]
			else:
				pos="complement("+s[3]+".."+s[4]+")"

			out=pos+"\tproduct\t"+prod

			if s[0] in rRNAdata:
				rRNAdata[s[0]]+="|"+out

			else:
				rRNAdata[s[0]]=out

	return rRNAdata

def split_gap(cds,gaps):	
	#single CDS		#cds:a..b; gaps:a..b|c..d|... → a..>egap<f..b
	if not "," in cds:	
		gs=[]		#[a,e>gap,<f,e>gap,<f,...,b]
		a,b=cds.split("..")
		aa=a.replace("<","");bb=b.replace(">","")

		aa=int(aa);bb=int(bb)
		gs.append(a)
		for gap in gaps.split("|"):
			e,f=gap.split("..")
		e=int(e);f=int(f)
		if e>bb or f<aa:
			pass
	#	elif e<aa:
	#		gs.remove(a);gs.append("<"+str(f+1))
	#	elif bb>f:
	#		b=">"+str(e-1)
		else:		
			e=">"+str(e-1)+"gap";f="<"+str(f+1)
			gs.append(e);gs.append(f)
		gs.append(b)
		if len(gs)==2:		#no gap
			return cds
		else:
			# remove cds length<100
			tmp=[]
			for i in range(1,len(gs),2):		
				i=i-1
				m=gs[i].replace("<","")
				n=gs[i+1].replace(">","").replace("gap","")
				if not int(n)-int(m)<100:
					tmp.append(gs[i]);tmp.append(gs[i+1])
			gs=tmp

			#[a,>egap,<f,>egap,<f,...,b] → [a,"..",<egap,",",f>,...,"..",b]
			c=["..",","]*((len(gs)-2)/2);c.append("..")
			cds2= [None]*(len(gs)+len(c))
			cds2[::2]=gs
			cds2[1::2]=c
			cds2="".join(str(x) for x in cds2)
	#multiple CDSs
	else:
		cds2=[]
		abs_list=cds.split(",")	#abs_list: [a..b,a..b,...]
		for ab in abs_list:
			gs=[]
			if len(ab.split(".."))==1:
				cds2.append(ab)
				continue
			a,b=ab.split("..")
			aa=a.replace("<","");bb=b.replace(">","")
			aa=int(aa);bb=int(bb)	
			gs.append(a)
			flg=False
			for gap in gaps.split("|"):
				e,f=gap.split("..")
				e=int(e);f=int(f)
				if e>bb or f<aa:
					pass
				elif e<aa:
					gs.remove(a);gs.append(f)
				elif bb<f:
					edge=">"+str(e-1)
					flg=True
				else:
					e=">"+str(e-1)+"gap";f="<"+str(f+1)
					gs.append(e);gs.append(f)
			if flg:
				gs.append(edge)
			else:
				gs.append(b)
			
			cds2+=gs

		# filter len<100
		gap=False
		for i in cds2:
			if i.endswith("gap"):
				gap=True
				ind=cds2.index(i)
		def calc_len(gs):
			length=0
			for i in range(2,len(gs)+1,2):
				i=i-1
				a=gs[i].replace("gap","").replace(">","")
				b=gs[i-1].replace("<","")
				a=int(a);b=int(b)
				length+=a-b
			return length
		if gap:
			cds21=cds2[:ind+1];cds22=cds2[ind+1:]
			if calc_len(cds21) < 100:
				cds2=[i for i in cds2 if i not in cds21]
			if calc_len(cds22) < 100:
				cds2=[i for i in cds2 if i not in cds22]

		else:
			if calc_len(cds2) < 100:
				cds2="";return cds2
		
		# add ".." and ","
		c=["..",","]*((len(cds2)-2)/2);c.append("..")	
		d= [None]*(len(cds2)+len(c))
		if (len(cds2) % 2) != 0:
			c.insert(0,",");d.append(None)
		d[::2]=cds2;d[1::2]=c;d="".join(str(x) for x in d)
		cds2=d

	return cds2

def quality_check(cds):
	quality="high"
	cds=cds.replace(">","");cds=cds.replace("<","")
	cds=cds.replace("join(","");cds=cds.replace(")","")
	a=re.split(',|\.\.',cds)
	for i in range(1,len(a)-1,2):
		if not int(a[i+1])-int(a[i])>15:
			quality="low"
			return quality
	return quality

def main():

	sps="Aspergillus udagawae"
	mLen=200

	header()

	##
	gff=sys.argv[1]
	funcFile=sys.argv[2]
	genomeFasta=sys.argv[3]
	rRNAmmer=sys.argv[4]

	## parse genome sequence fasta
	record_dict = SeqIO.to_dict(SeqIO.parse(genomeFasta, "fasta"))

	## get tRNA
	tRNA_data=gettRNA(gff)

	## get rRNA
	rRNA_data=getrRNA(rRNAmmer)
 
	## gap
	gapData=getGap(genomeFasta)

	##
	f=open(gff)
	lines=[line.rstrip() for line in f if line.startswith("NODE")]
	f.close()
	li1=[list(g) for k,g in groupby(lines, key=lambda x: x.split("\t")[0])]
	protfun=getProt(funcFile)

	##sort li1
	li1.sort(key=lambda x: int(x[0].split("_")[1]))

	rRNANum=1
	for scaffold in li1:
		scf=scaffold[0].split("\t")[0]	#
		length=int(scaffold[0].split("\t")[4])-int(scaffold[0].split("\t")[3])+1

		if length < mLen:
			continue

		scfName="scaffold_"+scf.split("_")[1]
		print("%s\tsource\t1..%s\torganism\tAspergillus udagawae" % (scfName,length))
		print("""\t\t\tstrain\tIFM 46972
\t\t\tmol_type\tgenomic DNA
\t\t\ttype_material\ttype strain
\t\t\tsubmitter_seqid\t@@[entry]@@
\t\t\tff_definition\t@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@
\t\t\tculture_collection\tIFM<JPN>:46972""")

		## check Gap information
		if scf in gapData:
			s=gapData[scf].split("|")

			for item in s:
				print("\tassembly_gap\t"+item+"\testimated_length\tknown")
				print("\t\t\tgap_type\twithin scaffold")
				print("\t\t\tlinkage_evidence\tpaired-ends")

		## check tRNA information
		if scf in tRNA_data:
			s=tRNA_data[scf].split("|")

			for item in s:
				if not "Pseudo" in item:
					print(item)

		##
		if scf in rRNA_data:
			s=rRNA_data[scf].split("|")

			for item in s:
				print("\trRNA\t"+item)
				out="\t\t\tlocus_tag\tIFM46972_r"+str(rRNANum)
				print(out)
				rRNANum+=1

		th=0
		for j in scaffold:
			if j.split("\t")[2]=="CDS":
				th=scaffold.index(j)
				break
		j=int(th)
	        if not j:
			continue
		a=scaffold[j].split("\t")

		strand=a[6]
		phase=a[7]

		while j <=len(scaffold):
			i=j;se=[]	#start&end

			try:
				while scaffold[i].split("\t")[2]=="CDS":
					a=scaffold[i].split("\t")
					feature=a[2];attribs=a[8].split(";");strand=a[6]
					aud=re.split("=|:",attribs[0])[1]	#locus tag
					se.append([a[3],a[4]])
					i+=1
			except IndexError:
				break

			## check start and end codon
			seq=record_dict[scf].seq

			if len(se)>=1:
				se=codonCheck(se, seq, strand)

			## single CDS
			if len(se)==1:
				cds= "%s..%s" % (se[0][0],se[0][1])
				
				if scf in gapData:
					cds2=split_gap(cds,gaps)
				else:
					cds2=cds

			## multiple CDSs
			if len(se) > 1:
				m=[]
				for k in se:
					if k[0]==k[1]:
						a=k[0]
					else:
						a=k[0]+".."+k[1]
					m.append(a)
				if strand=="-":
					m.reverse()
				cds=",".join(m)
				cds2=""
				if scf in gapData:
					gaps=gapData[scf]
					
					cds2=split_gap(cds,gaps)
					if cds2.startswith(","):
						cds2=cds2.replace(",","",1)
					cds2="join("+cds2+")"
				else:
					cds2="join("+cds+")"
			
			product="hypothetical protein"
			if aud in protfun:
				product=protfun[aud]
				if product=="":
					product="hypothetical protein"
			
			#move to note
			if product.endswith(")"):
				pro,note,tmp=re.split('\(|\)',product)
				product=pro+"\n\t\t\tnote\t"+note

			if se:	
				if not "gap," in cds2:
					if not cds2:
						continue
					cds2=cds2.replace("gap","")
					if strand=="+":
						print("\tCDS\t"+cds2+"\tproduct\t"+product)
					elif strand=="-":
						print("\tCDS\t"+"complement("+cds2+")\tproduct\t"+product)
					if not phase=="0":
						print("\t\t\tcodon_start\t"+phase)
					print("\t\t\tlocus_tag\t"+aud.split("-")[0])
					quality=quality_check(cds2)
					if quality=="low":
						print("\t\t\tartificial_location\tlow-quality sequence region")

				else:	#split cds by gap
					cds21,cds22=cds2.split("gap,")	
					if "join" in cds21:
						if "," in cds21:
							cds21=cds21+")"
						else:
							cds21=cds21.split("join(")[1]
						if "," in cds22:
							cds22="join("+cds22
						else:
							cds22=cds22.split(")")[0]
					
						
					for c in [cds21,cds22]:
						if strand=="+":
							print("\tCDS\t"+c+"\tproduct\t"+product)	
						else:
							print("\tCDS\t"+"complement("+c+")\tproduct\t"+product)		
						print("\t\t\tlocus_tag\t"+aud.split("-")[0])
						
						quality=quality_check(c)
						if quality=="low":
							print("\t\t\tartificial_location\tlow-quality sequence region")
					
				j=i#;print se
			else:
				j+=1

if __name__ == "__main__":
    main()	
