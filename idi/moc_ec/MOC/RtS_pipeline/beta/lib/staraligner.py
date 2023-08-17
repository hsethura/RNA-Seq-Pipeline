class StarAligner:
    def __init__(self, confd, outdir, refpath, corenum = 1):
        self.confd = confd
        self.outdir = outdir
        self.refpath = refpath
        self.corenum = corenum

    def align_paired(sample_id, read1, read2):
        confd = self.confd
        star_path = confd.STAR
        outdir = self.outdir
        refpath = self.refpath
        ldelim = confd.ldelim
        corenum = self.corenum

        out_prefix = outdir + ldelim + sample_id

        star_cmd = star_path + " --genomeDir " + refpath + " --readFilesIn " +
            read1 + " " + read2 + " --outFileNamePrefix " +  out_prefix +
            " --runThreadN " + corenum 
        call(star_cmd.split())
        print("Call: " + star_cmd)
        sam_path = out_prefix + "Aligned.out.sam"
        return sam_path

    def align_single(sample_id, read):
        confd = self.confd
        star_path = confd.STAR
        outdir = self.outdir
        refpath = self.refpath
        ldelim = confd.ldelim
        corenum = self.corenum

        out_prefix = outdir + ldelim + sample_id

        star_cmd = star_path + " --genomeDir " + refpath + " --readFilesIn " +
            read + " --outFileNamePrefix " +  out_prefix +
            " --runThreadN " + corenum 
        call(star_cmd.split())
        print("Call: " + star_cmd)
        sam_path = out_prefix + "Aligned.out.sam"
        return sam_path

