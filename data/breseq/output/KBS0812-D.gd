#=GENOME_DIFF	1.0
#=CREATED	05:35:50 15 Aug 2019
#=PROGRAM	breseq 0.32.0 revision 6ff6de7d1b87
#=COMMAND	breseq -j 8 -p -o /N/dc2/projects/muri2/Task2/LTDE/data/breseq/KBS0812-D -r /N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/KBS0812/GCF_005937795_2_ASM593779v2_genomic.gbff /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-D_S70_R1_001_clean.fastq.gz /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-D_S70_R2_001_clean.fastq.gz
#=REFSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/KBS0812/GCF_005937795_2_ASM593779v2_genomic.gbff
#=READSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-D_S70_R1_001_clean.fastq.gz
#=READSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-D_S70_R2_001_clean.fastq.gz
#=CONVERTED-BASES	606881971
#=CONVERTED-READS	9055396
#=INPUT-BASES	608071092
#=INPUT-READS	9071350
#=MAPPED-BASES	581126248
#=MAPPED-READS	8677248
SNP	1	5	NZ_CP041757	494250	A	frequency=1.36022568e-01
INS	2	6	NZ_CP041757	1755447	A	frequency=1	insert_position=1	repeat_length=1	repeat_new_copies=8	repeat_ref_copies=7	repeat_seq=A
SNP	3	7	NZ_CP041757	3176541	A	frequency=7.04221249e-01
SNP	4	8	NZ_CP041757	3715199	T	frequency=7.54232407e-02
RA	5	.	NZ_CP041757	494250	0	G	A	bias_e_value=2226580	bias_p_value=0.517828	consensus_score=159.7	fisher_strand_p_value=0.197488	frequency=1.36022568e-01	ks_quality_p_value=1	major_base=G	major_cov=44/28	major_frequency=8.63977432e-01	minor_base=A	minor_cov=10/2	new_cov=10/2	polymorphism_frequency=1.36022568e-01	polymorphism_score=10.9	prediction=polymorphism	ref_cov=44/28	total_cov=57/39
RA	6	.	NZ_CP041757	1755440	1	.	A	bias_e_value=4299850	bias_p_value=1	consensus_score=182.6	fisher_strand_p_value=1	frequency=1	ks_quality_p_value=1	major_base=A	major_cov=34/20	major_frequency=9.64358330e-01	minor_base=.	minor_cov=1/1	new_cov=34/20	polymorphism_frequency=9.64358330e-01	polymorphism_reject=SCORE_CUTOFF,FREQUENCY_CUTOFF,VARIANT_STRAND_COVERAGE,INDEL_HOMOPOLYMER	polymorphism_score=-1.0	prediction=consensus	ref_cov=1/1	total_cov=35/21
RA	7	.	NZ_CP041757	3176541	0	G	A	bias_e_value=4299850	bias_p_value=1	consensus_score=126.9	fisher_strand_p_value=1	frequency=7.04221249e-01	ks_quality_p_value=1	major_base=A	major_cov=54/15	major_frequency=7.04221249e-01	minor_base=G	minor_cov=23/6	new_cov=54/15	polymorphism_frequency=7.04221249e-01	polymorphism_score=71.3	prediction=polymorphism	ref_cov=23/6	total_cov=77/21
RA	8	.	NZ_CP041757	3715199	0	G	T	bias_e_value=10648.3	bias_p_value=0.00247644	consensus_score=633.0	fisher_strand_p_value=0.000268519	frequency=7.54232407e-02	ks_quality_p_value=1	major_base=G	major_cov=83/108	major_frequency=9.24576759e-01	minor_base=T	minor_cov=16/2	new_cov=16/2	polymorphism_frequency=7.54232407e-02	polymorphism_score=7.2	prediction=polymorphism	ref_cov=83/108	total_cov=99/110
JC	9	.	NZ_CP041757	1	1	NZ_CP041757	4215631	-1	0	alignment_overlap=0	circular_chromosome=1	coverage_minus=48	coverage_plus=65	flanking_left=75	flanking_right=75	frequency=1	junction_possible_overlap_registers=56	key=NZ_CP041757__1__1__NZ_CP041757__4215631__-1__0____75__75__0__0	max_left=68	max_left_minus=68	max_left_plus=64	max_min_left=35	max_min_left_minus=35	max_min_left_plus=31	max_min_right=37	max_min_right_minus=29	max_min_right_plus=37	max_pos_hash_score=112	max_right=66	max_right_minus=66	max_right_plus=63	neg_log10_pos_hash_p_value=NT	new_junction_coverage=1.05	new_junction_read_count=116	polymorphism_frequency=1.00000000e+00	pos_hash_score=45	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=56	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=56	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=113
JC	10	.	NZ_CP041758	1	1	NZ_CP041758	84215	-1	0	alignment_overlap=0	circular_chromosome=1	coverage_minus=115	coverage_plus=176	flanking_left=75	flanking_right=75	frequency=1	junction_possible_overlap_registers=56	key=NZ_CP041758__1__1__NZ_CP041758__84215__-1__0____75__75__0__0	max_left=68	max_left_minus=67	max_left_plus=68	max_min_left=36	max_min_left_minus=36	max_min_left_plus=35	max_min_right=35	max_min_right_minus=33	max_min_right_plus=35	max_pos_hash_score=112	max_right=69	max_right_minus=69	max_right_plus=67	neg_log10_pos_hash_p_value=NT	new_junction_coverage=1.16	new_junction_read_count=294	polymorphism_frequency=1.00000000e+00	pos_hash_score=67	prediction=consensus	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=0.00	side_1_overlap=0	side_1_possible_overlap_registers=56	side_1_read_count=0	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=0.00	side_2_overlap=0	side_2_possible_overlap_registers=56	side_2_read_count=0	side_2_redundant=0	total_non_overlap_reads=291
UN	11	.	NZ_CP041757	9356	14374
UN	12	.	NZ_CP041757	14377	14377
UN	13	.	NZ_CP041757	29908	34623
UN	14	.	NZ_CP041757	34892	34895
UN	15	.	NZ_CP041757	90097	94895
UN	16	.	NZ_CP041757	95422	95479
UN	17	.	NZ_CP041757	95597	95668
UN	18	.	NZ_CP041757	95919	100775
UN	19	.	NZ_CP041757	100778	100778
UN	20	.	NZ_CP041757	160516	165235
UN	21	.	NZ_CP041757	165725	165776
UN	22	.	NZ_CP041757	165993	170847
UN	23	.	NZ_CP041757	170995	175869
UN	24	.	NZ_CP041757	384373	384380
UN	25	.	NZ_CP041757	384513	384547
UN	26	.	NZ_CP041757	384667	384726
UN	27	.	NZ_CP041757	385389	385405
UN	28	.	NZ_CP041757	385673	385673
UN	29	.	NZ_CP041757	385676	385749
UN	30	.	NZ_CP041757	385914	385950
UN	31	.	NZ_CP041757	395135	395146
UN	32	.	NZ_CP041757	395284	395313
UN	33	.	NZ_CP041757	395441	395483
UN	34	.	NZ_CP041757	396153	396185
UN	35	.	NZ_CP041757	396440	396513
UN	36	.	NZ_CP041757	396685	396705
UN	37	.	NZ_CP041757	556351	556372
UN	38	.	NZ_CP041757	556497	556522
UN	39	.	NZ_CP041757	556725	556725
UN	40	.	NZ_CP041757	556729	556918
UN	41	.	NZ_CP041757	635168	639782
UN	42	.	NZ_CP041757	946245	948272
UN	43	.	NZ_CP041757	948396	950864
UN	44	.	NZ_CP041757	951606	951644
UN	45	.	NZ_CP041757	1870751	1870813
UN	46	.	NZ_CP041757	1967320	1967352
UN	47	.	NZ_CP041757	1968127	1968258
UN	48	.	NZ_CP041757	1968262	1968262
UN	49	.	NZ_CP041757	1968385	1968386
UN	50	.	NZ_CP041757	1968388	1968402
UN	51	.	NZ_CP041757	1968523	1968557
UN	52	.	NZ_CP041757	1978162	1978190
UN	53	.	NZ_CP041757	1978921	1979254
UN	54	.	NZ_CP041757	1979567	1979585
UN	55	.	NZ_CP041757	1979924	1979924
UN	56	.	NZ_CP041757	1979935	1980036
UN	57	.	NZ_CP041757	1980451	1980468
UN	58	.	NZ_CP041757	1980589	1980653
UN	59	.	NZ_CP041757	1993537	1993617
UN	60	.	NZ_CP041757	1994347	1994755
UN	61	.	NZ_CP041757	1995372	1995472
UN	62	.	NZ_CP041757	1995863	1995893
UN	63	.	NZ_CP041757	1996018	1996018
UN	64	.	NZ_CP041757	1996020	1996076
UN	65	.	NZ_CP041757	2161247	2161307
UN	66	.	NZ_CP041757	2246379	2246573
UN	67	.	NZ_CP041757	2246775	2246800
UN	68	.	NZ_CP041757	2246925	2246941
UN	69	.	NZ_CP041757	2579128	2579128
UN	70	.	NZ_CP041757	3172169	3172207
UN	71	.	NZ_CP041757	3172880	3172881
UN	72	.	NZ_CP041757	3172883	3172932
UN	73	.	NZ_CP041757	3173459	3176524
UN	74	.	NZ_CP041757	3176721	3178084
UN	75	.	NZ_CP041757	3178220	3178266
UN	76	.	NZ_CP041757	3666097	3666119
UN	77	.	NZ_CP041757	3673115	3673137
UN	78	.	NZ_CP041757	3855900	3855987
UN	79	.	NZ_CP041757	3855989	3855989
UN	80	.	NZ_CP041757	3856150	3856245
