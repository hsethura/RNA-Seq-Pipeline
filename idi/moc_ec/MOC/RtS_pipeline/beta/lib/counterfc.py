import re
from subprocess import call
from miscutils import get_pct

class CounterFC:
    def __init__(self, cldict, sampd):
        self.cldict = cldict
        self.sampd = sampd

  
    def get_reverse_dir(self, strand_dir):
        reverse_dir = ''
        if strand_dir == 'reverse':
            reverse_dir = 'forward'
        elif strand_dir == 'forward':
            reverse_dir = 'reverse'
        else:
            raise ValueError('Illegal value of strand_dir: ' + strand_dir)
        return reverse_dir

        
    def get_s_val_sense_allseq(self):
        cldict = self.cldict
        strand_dir = cldict.strand_dir
        s_val = ''
        if strand_dir == 'reverse':
            s_val = '2'
        elif strand_dir == 'forward':
            s_val = '1'
        return s_val

    def get_s_val_antisense_allseq(self):
        cldict = self.cldict
        strand_dir = cldict.strand_dir
        s_val = ''
        if strand_dir == 'reverse':
            s_val = '1'
        elif strand_dir == 'forward':
            s_val = '2'
        return s_val

    def get_s_val(self, exp_dir):
        s_val = ''
        if exp_dir == 's':
            s_val = self.get_s_val_sense_allseq()
        elif exp_dir == 'as':
            s_val = self.get_s_val_antisense_allseq()
        return s_val


    def exe_featureCounts(self, sample_id, ref_acc, outsorted, outdir, exp_dir):

        cldict = self.cldict
        ldelim = cldict.ldelim
        Data_dir = cldict.Data_dir

        s_val = self.get_s_val(exp_dir)
        print("Using s_val: " + s_val + " for exp_dir: " + exp_dir)

        featureCounts_str = cldict.featureCounts
        countfile_str = outdir + ldelim + sample_id + "_" + ref_acc + "_" + exp_dir + ".counts"
        patho_gff = Data_dir + ldelim + ref_acc + "_ALL.gff"
        fc_command_s = featureCounts_str + " -s " + s_val + " -t feature -g full_tag -Q 0 -a " + patho_gff + " -o " + countfile_str + " " +  outsorted
        print("Starting featureCounts..")
        print("fc_command_s: " + fc_command_s)
        call(fc_command_s.split())
        return countfile_str
       
    @staticmethod 
    def combine_s_as(countfile_s_str, countfile_as_str, countfile_str):
        inf_s = open(countfile_s_str, "r")
        inf_as = open(countfile_as_str, "r")
        outf = open(countfile_str, "w")
        # ignore first line
        inf_s.next()
        #outf.write("# Combining the content of sense and antisense reads\n")
        for line in inf_s:
            outf.write(line)
        inf_s.close()
        inf_as.next()
        inf_as.next()
        for line2 in inf_as:
            line3 = "AS_" + line2
            outf.write(line3)
        inf_as.close()
        outf.close()

    def load_metrics_from_sense_summary(self, sense_summary_str, metrics_str):
        """ This function would open this file and would calculate the 
        total number of reads by summing up the values on the right hand side.
        """
        sense_summary = open(sense_summary_str, "r")
        metrics = open(metrics_str, "w")
        # Ignore the first line 
        sense_summary.next()
        total_count = 0
        for line in sense_summary:
            parts = line.split()
            lcount = int(parts[-1])
            total_count += lcount
        metrics.write("metrics_type\tcount\n")
        metrics.write("total_count\t" + str(total_count) + "\n")
        sense_summary.close()
        metrics.close()
        return total_count


    def load_metrics_from_counts(self, counts_str, metrics_str, total_count):
        counts_file = open(counts_str, "r")
        metrics_file = open(metrics_str, "a")

        # Define the variables
        cds_count = 0
        igr_count = 0
        ncrna_count = 0
        rrna_count = 0
        trna_count = 0
        miscrna_count =0
        as_cds_count = 0
        as_igr_count = 0
        as_ncrna_count = 0
        as_rrna_count = 0
        as_trna_count = 0
        as_miscrna_count = 0

        counts_file.next()
        for line in counts_file:
            parts = line.split()
            tag = parts[0]
            count_str = parts[-1]
            lcount = int(count_str)
            is_as = None
            ltag = None
            if tag.startswith("AS_"):
                ltag = re.sub("^AS_", "", tag)
                is_as = True
            else:
                ltag = tag
                is_as = False
            
            if ltag.startswith("CDS:"):
                if is_as:
                    as_cds_count += lcount
                else:
                    cds_count += lcount
            elif ltag.startswith("rRNA:"):
                if is_as:
                    as_rrna_count += lcount
                else:
                    rrna_count += lcount
            elif ltag.startswith("tRNA:"):
                if is_as:
                    as_trna_count += lcount
                else:
                    trna_count += lcount
            elif ltag.startswith("ncRNA:"):
                if is_as:
                    as_ncrna_count += lcount
                else:
                    ncrna_count += lcount
            elif ltag.startswith("IGR:"):
                if is_as:
                    as_igr_count += lcount
                else:
                    igr_count += lcount
            elif ltag.startswith("misc_RNA:"):
                if is_as:
                    as_miscrna_count += lcount
                else:
                    miscrna_count += lcount
            else:
                raise ValueError('Unknown tag: ' + tag)

        sense_frag = cds_count + igr_count + rrna_count + \
            ncrna_count + trna_count + miscrna_count
        antisense_frag = as_cds_count + as_igr_count + as_rrna_count + \
            as_ncrna_count + as_trna_count + as_miscrna_count
        aligned_frag = sense_frag + antisense_frag
        unmapped_count = total_count - aligned_frag

        metrics_file.write("aligned_frag\t%d\n" % aligned_frag)
        aligned_frag_pct_of_total = get_pct(aligned_frag, total_count)
        metrics_file.write("aligned_frag_pct_of_total\t%.4f\n" % aligned_frag_pct_of_total)

        metrics_file.write("unmapped_count\t%d\n" % unmapped_count)
        unmapped_count_pct_of_total = get_pct(unmapped_count, total_count)
        metrics_file.write("unmapped_count_pct_of_total\t%.4f\n" % unmapped_count_pct_of_total)

        metrics_file.write("sense_frag\t%d\n" % sense_frag)
        sense_frag_pct_of_aligned = get_pct(sense_frag, aligned_frag)
        metrics_file.write("sense_frag_pct_of_aligned\t%.4f\n" % sense_frag_pct_of_aligned)

        self.print_feature_type("CDS", cds_count, as_cds_count, aligned_frag, metrics_file) 
        self.print_feature_type("IGR", igr_count, as_igr_count, aligned_frag, metrics_file) 
        self.print_feature_type("rRNA", rrna_count, as_rrna_count, aligned_frag, metrics_file) 
        self.print_feature_type("ncRNA", ncrna_count, as_ncrna_count, aligned_frag, metrics_file) 
        self.print_feature_type("tRNA", trna_count, as_trna_count, aligned_frag, metrics_file) 
        self.print_feature_type("miscRNA", miscrna_count, as_miscrna_count, aligned_frag, metrics_file) 
        counts_file.close()
        metrics_file.close()

    def print_feature_type(self, feature_name, count_s, count_as, aligned_frag, outfile):
        count_total = count_s + count_as
        feature_pct_of_aligned = get_pct(count_total, aligned_frag)
        feature_sense_pct = get_pct(count_s, count_total)
        outfile.write("%s_pct_of_aligned\t%.4f\n" % (feature_name, feature_pct_of_aligned))
        outfile.write("%s_sense_pct\t%.4f\n" % (feature_name, feature_sense_pct))

    def get_pct(numer, denom):
        if denom == 0:
            return 0.0
        else:
            lval = (numer * 100.0) / denom
            return lval

    def get_fc_metrics(self, sample_id, ref_acc, outdir):
        # This function would look into two files, once sense metric file
        # and another non-sense metric file and would create a final metric 
        # file similar to the one of the host side. Once that part is done
        # we can combine the metric files to make the unified metric file.

        cldict = self.cldict
        ldelim = cldict.ldelim
       
        sense_summary_str = outdir + ldelim + sample_id + "_" + ref_acc + \
            "_s.counts.summary" 
        counts_str = outdir + ldelim + sample_id + "_" + ref_acc + ".counts"
        metrics_str = outdir + ldelim + sample_id + "_" + ref_acc + ".metrics"

        print("Building metrics for  " + sample_id)
        print("sense_summary_str: " + sense_summary_str)
        print("counts_str: " + counts_str)
        print("metrics_str: " + metrics_str)

        total_count = self.load_metrics_from_sense_summary(sense_summary_str, 
            metrics_str)
        self.load_metrics_from_counts(counts_str, metrics_str, total_count)          
        
