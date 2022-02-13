import re, sys
def get_reverse_complementary(fa_seq):
	dic_gatc_com = {"G":"C","A":"T","T":"A","C":"G","g":"c","a":"t","t":"a","c":"g"}
	fa_seq_com_list = []
	for nuc in fa_seq:
		if nuc in dic_gatc_com.keys():
			fa_seq_com_list.append(dic_gatc_com[nuc])
		else:
			fa_seq_com_list.append(nuc)
	fa_seq_com_list.reverse()
	fa_seq_com = "".join(fa_seq_com_list)
	return fa_seq_com

def get_genome_fa(genome_fa_file):
	with open(genome_fa_file,"r") as genome_fa:
		dic_chr_seq = {}
		chrom = ""
		seq_list = [""]
		for line in genome_fa:
			if line.startswith(">"):
				dic_chr_seq[chrom] = "".join(seq_list)
				seq_list = []
				# chrom = line.strip()[1:]
				chrom = re.search("^>(\w+)",line.strip()).group(1)
			else:
				seq_list.append(line.strip())
		dic_chr_seq[chrom] = "".join(seq_list)
	genome_fa.close()
	return dic_chr_seq
def line_to_entry(line):
    f = line.rstrip().split("\t")
    gff_fields = f[:8]
    preattributes = re.split('\s*;\s*',f[8])
    attributes = {}
    for attribute in preattributes:
        m = re.search('(\S+)\s*["\']([^\'"]+)["\']',attribute)
        if m:
            attributes[m.group(1)] = m.group(2)
    if len(attributes.keys()) > 0:
        entry = {}
        entry['gff'] = gff_fields
        entry['attributes'] = attributes
        return entry
    return None
class GTFFile:
    def __init__(self,filename,genome_fa_file):
        self.genes = {}
        self.transcripts = {}
        with open(filename) as inf:
            for line in inf:
                if re.match('^#',line): continue
                f = line.rstrip().split("\t")
                if f[2] != 'exon': continue
                e = line_to_entry(line)
                if not e: continue
                if 'gene_id' not in e['attributes']:
                    sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n"+str(e)+"\n")
                    continue
                if 'transcript_id' not in e['attributes']:
                    sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n")
                    continue
                gene_id = e['attributes']['gene_id']
                transcript_id = e['attributes']['transcript_id']
                if gene_id not in self.genes:
                    self.genes[gene_id] = []
                self.genes[gene_id].append(e)
                if transcript_id not in self.transcripts:
                    self.transcripts[transcript_id] = []
                self.transcripts[transcript_id].append(e)
        self.dic_chr_seq = get_genome_fa(genome_fa_file)
        return
    def write_genepred(self,filehandle):
        for transcript in self.transcripts:
            #filehandle.write(str(transcript)+"\n")
            tlist = self.transcripts[transcript]
            elist = {}
            gene = '.'
            strand = '.'
            chrom = '.'
            for t in tlist:
                if not t['gff'][2].lower() == 'exon':
                    continue
                start = int(t['gff'][3])-1
                end = int(t['gff'][4])
                if start > end:
                    sys.stderr.write("ERROR start bigger than end\n")
                    sys.exit()
                elist[start] = end
                gene = t['attributes']['gene_id']
                strand = t['gff'][6]
                chrom = t['gff'][0]
            #print strand
            #print transcript
            starts = sorted(elist,key=lambda k: elist[k])
            first = starts[0]
            last = elist[starts[len(starts)-1]]
            ostring = gene + "\t" + transcript + "\t" + chrom + "\t" + strand + "\t" + str(first) + "\t" \
                            + str(last) + "\t" + str(first) + "\t" + str(last) + "\t" \
                            + str(len(starts)) + "\t" \
                            + ",".join([str(x) for x in starts]) + ",\t" \
                            + ",".join([str(elist[x]) for x in starts])+","
            filehandle.write(ostring+"\n")
    def write_fa(self,filehandle):
        dic_iso_fa = {}
        for transcript in self.transcripts:
            #filehandle.write(str(transcript)+"\n")
            tlist = self.transcripts[transcript]
            elist = {}
            gene = '.'
            strand = '.'
            chrom = '.'
            for t in tlist:
                if not t['gff'][2].lower() == 'exon':
                    continue
                start = int(t['gff'][3])-1
                end = int(t['gff'][4])
                if start > end:
                    sys.stderr.write("ERROR start bigger than end\n")
                    sys.exit()
                elist[start] = end
                gene = t['attributes']['gene_id']
                strand = t['gff'][6]
                chrom = t['gff'][0]
            if chrom not in self.dic_chr_seq:
                continue
            #print strand
            #print transcript
            starts = sorted(elist,key=lambda k: elist[k])
            first = starts[0]
            last = elist[starts[len(starts)-1]]
            seq_list = []
            for start in starts:
                end = int(elist[start])
                seq_list.append(self.dic_chr_seq[chrom][start:end])
            
            seq = "".join(seq_list)
            tss,tts,extend_both_end = first,last,0
            if strand == "+":
                if (int(tss)-int(extend_both_end)) > 0:
                    seq = self.dic_chr_seq[chrom][(int(tss)-int(extend_both_end)):int(tss)] + seq
                else:
                    seq = self.dic_chr_seq[chrom][0:int(tss)] + seq
                if (int(tts)+int(extend_both_end)) < (len(self.dic_chr_seq[chrom])-1):
                    seq = seq + self.dic_chr_seq[chrom][int(tts):(int(tts)+int(extend_both_end))]
                else:
                    seq = seq + self.dic_chr_seq[chrom][int(tts):(len(self.dic_chr_seq[chrom])-1)]
            else:
                if (int(tss)-int(extend_both_end)) > 0:
                    seq = self.dic_chr_seq[chrom][(int(tss)-int(extend_both_end)):int(tss)] + seq
                else:
                    seq = self.dic_chr_seq[chrom][0:int(tss)] + seq
                if (int(tts)+int(extend_both_end)) < (len(self.dic_chr_seq[chrom])-1):
                    seq = seq + self.dic_chr_seq[chrom][int(tts):(int(tts)+int(extend_both_end))]
                else:
                    seq = seq + self.dic_chr_seq[chrom][int(tts):(len(self.dic_chr_seq[chrom])-1)]
                seq = get_reverse_complementary(seq)
            dic_iso_fa[transcript] = seq
        with open(filehandle,'w') as f:
            for iso in dic_iso_fa:
                f.write(">" + iso + '\n')
                f.write(dic_iso_fa[iso]+ '\n')