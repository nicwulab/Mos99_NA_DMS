# only variable needed to change


PROJECT_PATH='/Users/yiquan/PycharmProjects/Mos99'
FASTA = PROJECT_PATH + '/result/{SAMPLENAME}_tagfree.fa'
SAMPLENAMES, = glob_wildcards(FASTA)
REF=PROJECT_PATH + '/ref/Mos99_NA.fa'
#print(SAMPLENAMES)
FAS=PROJECT_PATH + '/result/{SAMPLENAME}_tagfree.fa'
RESULT_PATH = PROJECT_PATH + '/result/{SAMPLENAME}'
RM_PRIMER_FAS = RESULT_PATH + '/{SAMPLENAME}_rm_primer.fa'
UNTRIM_FAS = RESULT_PATH + '/{SAMPLENAME}_untrim.fa'
RM_PRIMER_FQ = RESULT_PATH + '/{SAMPLENAME}_rm_primer.fq'
ASSEMBLED_FQ = RESULT_PATH + '/{SAMPLENAME}_assembled.fq'
TABLE = RESULT_PATH + '/{SAMPLENAME}_count.tsv'

rule all:
    input:
        expand(TABLE, SAMPLENAME=SAMPLENAMES),
        expand(RM_PRIMER_FAS, SAMPLENAME = SAMPLENAMES),
        expand(RM_PRIMER_FQ, SAMPLENAME = SAMPLENAMES),
        expand(ASSEMBLED_FQ, SAMPLENAME = SAMPLENAMES)
# remove primer(only for paired-primer are removed) &filter out small reads(-m 100:100(R1:R2))
rule rm_primer:
    input:
        FAS
    params:
        rename=lambda wc: "'{id}_{adapter_name}'"
    output:
        trim_o=RM_PRIMER_FAS,
        untrim_o = UNTRIM_FAS
    shell:
        '''cutadapt --pair-adapters --interleaved '''\
        '''-g "Amp1=AAGGAAATATGCCCCAAACTA;e=0.2" -g "Amp2=ACACTAAACAACGGGCATTCA;e=0.2" -g "Amp3=GCTAGCTTCATTTACAATGGG;e=0.2" -g "Amp4=CCATTGTCAGGAAGTGCTCAG;e=0.2" -g "Amp5=AGCTCCAGCAGTAGCCATTGC;e=0.2" -g "Amp6=CAAGTCATAGTTGACAGAGGT;e=0.2" '''\
        '''-G "Amp1=CCTATCATGTACTGTGTCATT;e=0.2" -G "Amp2=ACCAATACTATCTACAAGCCT;e=0.2" -G "Amp3=ACAGGAGCATTCCTCGACATG;e=0.2" -G "Amp4=TTCCTCATTGTTAGGATCCAA;e=0.2" -G "Amp5=ACCAGAATAACCGGACCTATT;e=0.2" -G "Amp6=GAACAAATTATATAGGCATGAG;e=0.2" '''\
        '''--rename={params} '''\
        '''--untrimmed-output {output.untrim_o} '''
        '''-O 10 -m 50:50 -o {output.trim_o} {input}'''

rule fas2fq:
    input:
        RM_PRIMER_FAS
    output:
        RM_PRIMER_FQ
    shell:
        "seqtk seq -F '#' {input} > {output}"
#merge paired-end reads
rule flash:
    input:
        RM_PRIMER_FQ
    output:
        ASSEMBLED_FQ
    shell:
        "flash -m 30 -M 70 -I {input} -c > {output} "

rule fq2count:
    input:
        ASSEMBLED_FQ
    params:
        REF_FA=REF
    output:
        TABLE
    shell:
        'python fastq2count.py {input} {params} {output}'
