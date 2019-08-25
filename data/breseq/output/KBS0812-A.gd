#=GENOME_DIFF	1.0
#=CREATED	09:05:31 25 Jun 2019
#=PROGRAM	breseq 0.32.0 revision 6ff6de7d1b87
#=COMMAND	breseq -j 8 -p -o /N/dc2/projects/muri2/Task2/LTDE/data/breseq/KBS0812-A -r /N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/KBS0812/GCF_005937795.1_ASM593779v1_genomic.gbff /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-A_S49_R1_001_clean.fastq.gz /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-A_S49_R2_001_clean.fastq.gz
#=REFSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/KBS0812/GCF_005937795.1_ASM593779v1_genomic.gbff
#=READSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-A_S49_R1_001_clean.fastq.gz
#=READSEQ	/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-KBS0812-A_S49_R2_001_clean.fastq.gz
#=CONVERTED-BASES	574157740
#=CONVERTED-READS	8574661
#=INPUT-BASES	575382758
#=INPUT-READS	8591090
#=MAPPED-BASES	556404026
#=MAPPED-READS	8314976
SNP	1	9	NZ_VCBS01000001	44919	T	frequency=5.47943115e-02
INS	2	10	NZ_VCBS01000001	601081	A	frequency=1	insert_position=1	repeat_length=1	repeat_new_copies=8	repeat_ref_copies=7	repeat_seq=A
SNP	3	11	NZ_VCBS01000001	647524	G	frequency=6.12578392e-02
SNP	4	12	NZ_VCBS01000001	2562738	T	frequency=2.18108654e-01
SNP	5	13	NZ_VCBS01000002	393498	A	frequency=8.95013809e-02
SNP	6	14	NZ_VCBS01000002	393522	A	frequency=8.14738274e-02
SNP	7	15	NZ_VCBS01000002	487200	A	frequency=6.82506561e-02
SNP	8	16	NZ_VCBS01000002	950095	A	frequency=5.87739944e-02
RA	9	.	NZ_VCBS01000001	44919	0	G	T	bias_e_value=4286250	bias_p_value=0.998974	consensus_score=283.0	fisher_strand_p_value=1	frequency=5.47943115e-02	ks_quality_p_value=0.955041	major_base=G	major_cov=54/32	major_frequency=9.45205688e-01	minor_base=T	minor_cov=3/2	new_cov=3/2	polymorphism_frequency=5.47943115e-02	polymorphism_score=4.8	prediction=polymorphism	ref_cov=54/32	total_cov=57/34
RA	10	.	NZ_VCBS01000001	601074	1	.	A	consensus_score=278.0	frequency=1	major_base=A	major_cov=58/24	major_frequency=1.00000000e+00	minor_base=T	minor_cov=0/1	new_cov=58/24	polymorphism_frequency=1.00000000e+00	polymorphism_reject=SCORE_CUTOFF,FREQUENCY_CUTOFF,VARIANT_STRAND_COVERAGE,SURROUNDING_HOMOPOLYMER	polymorphism_score=-6.6	prediction=consensus	ref_cov=1/0	total_cov=59/25
RA	11	.	NZ_VCBS01000001	647524	0	A	G	bias_e_value=1967690	bias_p_value=0.4586	consensus_score=285.6	fisher_strand_p_value=0.162954	frequency=6.12578392e-02	ks_quality_p_value=1	major_base=A	major_cov=21/69	major_frequency=9.38742161e-01	minor_base=G	minor_cov=3/3	new_cov=3/3	polymorphism_frequency=6.12578392e-02	polymorphism_score=2.7	prediction=polymorphism	ref_cov=21/69	total_cov=25/72
RA	12	.	NZ_VCBS01000001	2562738	0	G	T	bias_e_value=0.000527925	bias_p_value=1.23041e-10	consensus_score=303.2	fisher_strand_p_value=4.53711e-12	frequency=2.18108654e-01	ks_quality_p_value=1	major_base=G	major_cov=89/34	major_frequency=7.81891346e-01	minor_base=T	minor_cov=2/30	new_cov=2/30	polymorphism_frequency=2.18108654e-01	polymorphism_score=10.1	prediction=polymorphism	ref_cov=89/34	total_cov=96/93
RA	13	.	NZ_VCBS01000002	393498	0	C	A	bias_e_value=4290650	bias_p_value=1	consensus_score=537.9	fisher_strand_p_value=1	frequency=8.95013809e-02	ks_quality_p_value=1	major_base=C	major_cov=119/62	major_frequency=9.10498619e-01	minor_base=A	minor_cov=17/8	new_cov=17/8	polymorphism_frequency=8.95013809e-02	polymorphism_score=6.2	prediction=polymorphism	ref_cov=119/62	total_cov=136/71
RA	14	.	NZ_VCBS01000002	393522	0	T	A	bias_e_value=68578.7	bias_p_value=0.0159833	consensus_score=377.5	fisher_strand_p_value=0.00225253	frequency=8.14738274e-02	ks_quality_p_value=1	major_base=T	major_cov=72/68	major_frequency=9.18526173e-01	minor_base=A	minor_cov=16/2	new_cov=16/2	polymorphism_frequency=8.14738274e-02	polymorphism_score=2.1	prediction=polymorphism	ref_cov=72/68	total_cov=88/70
RA	15	.	NZ_VCBS01000002	487200	0	C	A	bias_e_value=111629	bias_p_value=0.0260169	consensus_score=397.9	fisher_strand_p_value=0.00398751	frequency=6.82506561e-02	ks_quality_p_value=1	major_base=C	major_cov=46/75	major_frequency=9.31749344e-01	minor_base=A	minor_cov=10/2	new_cov=10/2	polymorphism_frequency=6.82506561e-02	polymorphism_score=5.6	prediction=polymorphism	ref_cov=46/75	total_cov=56/77
RA	16	.	NZ_VCBS01000002	950095	0	C	A	bias_e_value=5041.37	bias_p_value=0.00117496	consensus_score=476.2	fisher_strand_p_value=0.000116859	frequency=5.87739944e-02	ks_quality_p_value=1	major_base=C	major_cov=111/32	major_frequency=9.41226006e-01	minor_base=A	minor_cov=2/9	new_cov=2/9	polymorphism_frequency=5.87739944e-02	polymorphism_score=4.6	prediction=polymorphism	ref_cov=111/32	total_cov=114/41
MC	17	.	NZ_VCBS01000004	1	60	0	58	left_inside_cov=0	left_outside_cov=NA	right_inside_cov=52	right_outside_cov=61
JC	18	.	NZ_VCBS01000001	1984037	1	NZ_VCBS01000001	1984046	1	0	alignment_overlap=2	coverage_minus=0	coverage_plus=3	flanking_left=75	flanking_right=75	frequency=2.10204472e-02	junction_possible_overlap_registers=53	key=NZ_VCBS01000001__1984037__1__NZ_VCBS01000001__1984044__1__2____75__75__0__0	max_left=53	max_left_minus=0	max_left_plus=53	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=36	max_min_right_minus=0	max_min_right_plus=36	max_pos_hash_score=108	max_right=36	max_right_minus=0	max_right_plus=36	neg_log10_pos_hash_p_value=NT	new_junction_coverage=0.03	new_junction_read_count=3	polymorphism_frequency=2.10204472e-02	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=1.48	side_1_overlap=2	side_1_possible_overlap_registers=55	side_1_read_count=153	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=1.32	side_2_overlap=0	side_2_possible_overlap_registers=53	side_2_read_count=132	side_2_redundant=0	total_non_overlap_reads=3
JC	19	.	NZ_VCBS01000001	3076736	1	NZ_VCBS01000001	3076745	1	0	alignment_overlap=4	coverage_minus=0	coverage_plus=4	flanking_left=75	flanking_right=75	frequency=3.80293863e-02	junction_possible_overlap_registers=51	key=NZ_VCBS01000001__3076736__1__NZ_VCBS01000001__3076741__1__4____75__75__0__0	max_left=45	max_left_minus=0	max_left_plus=45	max_min_left=0	max_min_left_minus=0	max_min_left_plus=0	max_min_right=17	max_min_right_minus=0	max_min_right_plus=17	max_pos_hash_score=104	max_right=17	max_right_minus=0	max_right_plus=17	neg_log10_pos_hash_p_value=NT	new_junction_coverage=0.04	new_junction_read_count=4	polymorphism_frequency=3.80293863e-02	pos_hash_score=3	prediction=polymorphism	reject=FREQUENCY_CUTOFF	side_1_annotate_key=gene	side_1_continuation=0	side_1_coverage=1.01	side_1_overlap=4	side_1_possible_overlap_registers=55	side_1_read_count=105	side_1_redundant=0	side_2_annotate_key=gene	side_2_continuation=0	side_2_coverage=1.09	side_2_overlap=0	side_2_possible_overlap_registers=51	side_2_read_count=105	side_2_redundant=0	total_non_overlap_reads=4
UN	20	.	NZ_VCBS01000001	1	8
UN	21	.	NZ_VCBS01000001	716382	716450
UN	22	.	NZ_VCBS01000001	812962	812985
UN	23	.	NZ_VCBS01000001	813766	813888
UN	24	.	NZ_VCBS01000001	814027	814036
UN	25	.	NZ_VCBS01000001	814157	814182
UN	26	.	NZ_VCBS01000001	823803	823824
UN	27	.	NZ_VCBS01000001	824562	824669
UN	28	.	NZ_VCBS01000001	824757	824881
UN	29	.	NZ_VCBS01000001	825577	825666
UN	30	.	NZ_VCBS01000001	826085	826102
UN	31	.	NZ_VCBS01000001	826227	826292
UN	32	.	NZ_VCBS01000001	839175	839251
UN	33	.	NZ_VCBS01000001	839967	840048
UN	34	.	NZ_VCBS01000001	840184	840390
UN	35	.	NZ_VCBS01000001	841010	841091
UN	36	.	NZ_VCBS01000001	841506	841507
UN	37	.	NZ_VCBS01000001	841510	841529
UN	38	.	NZ_VCBS01000001	841654	841710
UN	39	.	NZ_VCBS01000001	1006881	1006943
UN	40	.	NZ_VCBS01000001	1092009	1092009
UN	41	.	NZ_VCBS01000001	1092011	1092214
UN	42	.	NZ_VCBS01000001	1092407	1092407
UN	43	.	NZ_VCBS01000001	1092409	1092434
UN	44	.	NZ_VCBS01000001	1092558	1092573
UN	45	.	NZ_VCBS01000001	2017803	2017841
UN	46	.	NZ_VCBS01000001	2018504	2018566
UN	47	.	NZ_VCBS01000001	2019093	2019191
UN	48	.	NZ_VCBS01000001	2019320	2019320
UN	49	.	NZ_VCBS01000001	2019330	2021808
UN	50	.	NZ_VCBS01000001	2021881	2022159
UN	51	.	NZ_VCBS01000001	2022356	2023718
UN	52	.	NZ_VCBS01000001	2023854	2023883
UN	53	.	NZ_VCBS01000001	2023886	2023886
UN	54	.	NZ_VCBS01000001	2511724	2511753
UN	55	.	NZ_VCBS01000001	2518751	2518772
UN	56	.	NZ_VCBS01000001	2701534	2701619
UN	57	.	NZ_VCBS01000001	2701623	2701623
UN	58	.	NZ_VCBS01000001	2701782	2701877
UN	59	.	NZ_VCBS01000001	3070632	3072227
UN	60	.	NZ_VCBS01000001	3072552	3072822
UN	61	.	NZ_VCBS01000001	3072900	3075615
UN	62	.	NZ_VCBS01000001	3091173	3091211
UN	63	.	NZ_VCBS01000002	1	17
UN	64	.	NZ_VCBS01000002	59761	64480
UN	65	.	NZ_VCBS01000002	64959	65018
UN	66	.	NZ_VCBS01000002	65233	70089
UN	67	.	NZ_VCBS01000002	70235	75107
UN	68	.	NZ_VCBS01000002	283613	283622
UN	69	.	NZ_VCBS01000002	283761	283788
UN	70	.	NZ_VCBS01000002	283915	283961
UN	71	.	NZ_VCBS01000002	284636	284647
UN	72	.	NZ_VCBS01000002	284649	284650
UN	73	.	NZ_VCBS01000002	284918	284918
UN	74	.	NZ_VCBS01000002	284921	284996
UN	75	.	NZ_VCBS01000002	285171	285183
UN	76	.	NZ_VCBS01000002	294383	294386
UN	77	.	NZ_VCBS01000002	294526	294546
UN	78	.	NZ_VCBS01000002	294673	294723
UN	79	.	NZ_VCBS01000002	295400	295414
UN	80	.	NZ_VCBS01000002	295686	295749
UN	81	.	NZ_VCBS01000002	295936	295947
UN	82	.	NZ_VCBS01000002	455600	455614
UN	83	.	NZ_VCBS01000002	455739	455762
UN	84	.	NZ_VCBS01000002	455966	456162
UN	85	.	NZ_VCBS01000002	534414	539030
UN	86	.	NZ_VCBS01000002	845490	847520
UN	87	.	NZ_VCBS01000002	847597	850115
UN	88	.	NZ_VCBS01000002	850854	850897
UN	89	.	NZ_VCBS01000002	1053683	1053691
UN	90	.	NZ_VCBS01000003	1	10
UN	91	.	NZ_VCBS01000003	15	15
UN	92	.	NZ_VCBS01000003	18	18
UN	93	.	NZ_VCBS01000003	84275	84292
UN	94	.	NZ_VCBS01000004	1	13
UN	95	.	NZ_VCBS01000004	55215	56805
UN	96	.	NZ_VCBS01000004	56875	60013
UN	97	.	NZ_VCBS01000004	60537	60591
UN	98	.	NZ_VCBS01000004	60724	60786
UN	99	.	NZ_VCBS01000004	61042	61177
UN	100	.	NZ_VCBS01000005	1	16
UN	101	.	NZ_VCBS01000005	277	282
