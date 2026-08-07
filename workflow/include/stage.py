from collections import OrderedDict


class Stage:
    def __init__(self, config, stage_name, logger, prev_stage=None):
        self.stage_name = stage_name
        self.prev_stage = prev_stage
        self.logger = logger
        self.parameters = OrderedDict()
        self.config = config
        self.logger.info(f"Initializing stage {self.stage_name}...")

        self.no_fasta_stage_set = set(config["no_fasta_stage_list"])
        self.assembly_initiating_stage_set = set(config["assembly_initiating_stage_list"])
        """
        if stage_name == "contig":
            for assembler in config["stage_coretools"][self.stage_name]:
                #print(parameters["tool_options"][assembler])
                if assembler == "hifiasm":
                    option_set_group_dict, option_set_group_assignment_dict = group_option_sets(parameters["tool_options"][assembler],
                                                                                                config["tool_specific_features"][assembler]["options_affecting_error_correction"],
                                                                                                tool=assembler)

                else:
                    option_set_group_assignment_dict = None
                    option_set_group_dict = None
                assembler_option_set_group_dict[assembler] = option_set_group_dict
        """
        if stage_name not in self.no_fasta_stage_set:
            self.logger.info(TAB + f"Requested parameters:")
            for tool in config["stage_coretools"][self.stage_name]:
                for option_set in config["coretool_option_sets"][tool]:
                    parameters_label_list = []
                    if self.stage_name == "draft_qc": #in self.assembly_initiating_stage_set:
                        parameters_label_list.append("{0}_{1}".format(tool, option_set))
                    elif self.stage_name == "contig":
                        ploidy = self.detect_ploidy_for_contig_stage(tool, option_set)
                        parameters["tool_options"][tool][option_set]["assembly_ploidy"] = ploidy
                        parameters_label_list.append("{0}_{1}@p{2}".format(tool, option_set, ploidy))
                    elif self.stage_name == "dedup":
                        if tool in ["hapsolo", "combo_purge"]:
                            for busco_lineage in config["tool_manually_adjusted_features"]["hapsolo"]["busco_lineage_list"]:
                                for prev_parameters in stage_dict[self.prev_stage].parameters:
                                    parameters_label_list.append("{0}..{1}_{2}@{3}".format(prev_parameters, tool, option_set, busco_lineage)),
                        else:
                            for prev_parameters in stage_dict[self.prev_stage].parameters:
                                parameters_label_list.append("{0}..{1}_{2}".format(prev_parameters, tool, option_set))

                    elif self.stage_name == "ref_scaffolding":
                        if "reference" not in config["data"]:
                            raise ValueError("ERROR!!! ref_scaffolding stage was requested, but no reference was found!!! Check the main config file and input files...")
                        for reference in config["data"]["reference"]["ref_dict"]:
                            for prev_parameters in stage_dict[self.prev_stage].parameters:
                                parameters_label_list.append("{0}..{1}_{2}@{3}".format(prev_parameters, tool, option_set, reference))

                    elif self.stage_name == "hic_scaffolding":
                        if tool == "threeddna": # threeddna works only in connection with juicer
                            for prev_parameters in stage_dict[self.prev_stage].parameters:
                                if "juicer" in prev_parameters:
                                    parameters_label_list.append("{0}..{1}_{2}".format(prev_parameters, tool, option_set)),
                        else:
                            for prev_parameters in stage_dict[self.prev_stage].parameters:
                                if "juicer" in prev_parameters: # usually juicer is used with threeddna only, it is possible to use it with other hic scaffolders, but for what?
                                    continue
                                parameters_label_list.append("{0}..{1}_{2}".format(prev_parameters, tool, option_set))

                    elif self.prev_stage not in self.no_fasta_stage_set:
                        for prev_parameters in stage_dict[self.prev_stage].parameters:
                            parameters_label_list.append("{0}..{1}_{2}".format(prev_parameters, tool, option_set))
                    else:
                        raise ValueError(f"ERROR!!! Stage {self.stage_name} is not an assembly initiating stage (i.e. does not create or get fasta from input), "
                                         f"and previous stage ({self.prev_stage} is not an assembly initiating stage as well. Check 'stage_list' in your main config file )")

                    for parameters_label in parameters_label_list:
                        self.logger.info(TAB * 2 + parameters_label)
                        self.parameters[parameters_label] = {}

                        self.parameters[parameters_label]["tool"] = tool
                        self.parameters[parameters_label]["option_set"] = deepcopy(parameters["tool_options"][tool][option_set])
                        if "main_datatypes" in self.parameters[parameters_label]["option_set"]:
                            self.parameters[parameters_label]["option_set"]["main_datatypes"] = set(self.parameters[parameters_label]["option_set"]["main_datatypes"]) & set(self.config["data"].keys())
                            self.logger.info(TAB * 3 + f"Main datatypes: {', '.join(self.parameters[parameters_label]['option_set']['main_datatypes'])}")
                        if self.stage_name in self.assembly_initiating_stage_set:
                            self.parameters[parameters_label]["prev_parameters"] = None

                            # set ploidy of the assembly
                            if self.stage_name == "draft_qc":
                                self.parameters[parameters_label]["option_set"]["assembly_ploidy"] = len(config["data"]["draft"]["haplotypes"])
                            elif self.stage_name == "contig":
                                if tool == "hifiasm":
                                    option_set_group_dict, option_set_group_assignment_dict = group_option_sets(parameters["tool_options"][tool],
                                                                                                                config["tool_specific_features"][tool]["options_affecting_error_correction"],
                                                                                                                tool=tool)
                                else:
                                    option_set_group_assignment_dict = None
                                    option_set_group_dict = None

                                assembler_option_set_group_dict[tool] = option_set_group_dict

                                self.parameters[parameters_label]["option_set_group"] = option_set_group_assignment_dict[option_set] if option_set_group_assignment_dict is not None else None

                            else:
                                raise ValueError(f"ERROR!!! Unknown assembly initiating stage {self.stage_name}! Something is wrong...")
                            self.parameters[parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, self.parameters[parameters_label]["option_set"]["assembly_ploidy"] + 1)] if self.parameters[parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]

                        else:
                            # inherit haplotypes from previous stage
                            self.parameters[parameters_label]["prev_parameters"] = "..".join(parameters_label.split("..")[:-1])
                            self.parameters[parameters_label]["haplotype_list"] = stage_dict[self.prev_stage].parameters[self.parameters[parameters_label]["prev_parameters"]]["haplotype_list"]

                            if len(self.parameters[parameters_label]["haplotype_list"]) == 1:
                                if ("use_phased_reads" in self.parameters[parameters_label]["option_set"]) and self.parameters[parameters_label]["option_set"]["use_phased_reads"]:
                                    self.parameters.pop(parameters_label)
                                    self.logger.warn(f"WARNING!!! Impossible to phase reads for {parameters_label} of {self.stage_name} as input draft assembly is haploid! Excluding it...")
                        if parameters_label in self.parameters:
                            self.logger.info(TAB * 3 + f"Assembly ploidy: {len(self.parameters[parameters_label]['haplotype_list'])}")
                            self.logger.info(TAB * 3 + f"Haplotype list: {', '.join(self.parameters[parameters_label]['haplotype_list'])}")

    def detect_ploidy_for_contig_stage(self, tool, option_set):
        found_dict = {}
        if tool in ["hifiasm", "verkko"]:
            for datatype in "hic", "parental":
                if parameters["tool_options"][tool][option_set][f"use_{datatype}"]:
                    #self.logger.info(TAB * 3 + f"Parameter set allows usage of phasing datatype {datatype}...")
                    if datatype in self.config["data"]:
                        found_dict[datatype] = True
                        #self.logger.info(TAB * 4 + f"Datatype {datatype} was found...")
                    else:
                        found_dict[datatype] = False
                        #self.logger.info(TAB * 4 + f"Datatype {datatype} was not found...")
                #else:
                #    #self.logger.info(TAB * 3 + f"Parameter set doesn't allow usage of phasing datatype {datatype}...")

            if sum(found_dict.values()) >= 1:
                #self.logger.info(TAB * 3 + f"Phasing data allowed to use was found...")
                if ("assembly_ploidy" in parameters["tool_options"][tool][option_set]) and parameters["tool_options"][tool][option_set]["assembly_ploidy"]:
                    return parameters["tool_options"][tool][option_set]["assembly_ploidy"]
                    #self.logger.info(TAB * 3 + f"Parameter set enforces assembly ploidy {parameters["tool_options"][tool][option_set]['assembly_ploidy']}...")
                else:
                    #self.logger.info(TAB * 3 + f"Parameter doesn't enforce assembly ploidy. Setting it to {self.config['ploidy']} according to the main config file... ")
                    return self.config['ploidy']
            else:
                #self.logger.info(TAB * 3 + f"Phasing data allowed to use was not found...")
                #self.logger.info(TAB * 3 + f"Enforcing assembly ploidy 1...")
                return 1
        elif tool == "flye":
            return 1
        elif tool == "nextdenovo":
            return 1
        else:
            raise ValueError(f"ERROR!!! Impossible to set ploidy for tool {tool}. Check if it is a typo, or add code related to this tool to detect_ploidy_for_contig_stage method of Stage class from stages.py.")

    def request_files(self):
        results_list = []
        if self.stage_name not in self.no_fasta_stage_set | {"hic_alignment",}:
            results_list += self.request_common_fasta_stage_files()

        if self.stage_name == "raw_read_qc":
            results_list += self.request_qc_files(stage="raw")
        if self.stage_name == "raw_kmer_qc":
            results_list += self.request_kmer_qc_files(stage="raw")
        if self.stage_name == "filter_reads":
            results_list += self.request_filter_reads_files()
        if self.stage_name == "filtered_read_qc":
            results_list += self.request_qc_files(stage="filtered")
        if self.stage_name == "read_contamination_scan":
            results_list += self.request_read_contamination_scan_files()
        if self.stage_name == "mtdna":
            results_list += self.request_mtdna_files()
        if self.stage_name == "kmer_qc":
            results_list += self.request_kmer_qc_files(stage="final")
        if self.stage_name == "contig":
            results_list += self.request_contig_files()
        if self.stage_name == "hic_alignment":
            results_list += self.request_hic_alignment_files()
        if self.stage_name == "hic_scaffolding":
            results_list += self.request_hic_scaffolding_files()
        if self.stage_name == "read_phasing":
            results_list += self.request_read_phasing_files()
        if self.stage_name == "dedup":
            results_list += self.request_dedup_files()
        if self.stage_name == "ref_scaffolding":
            results_list += self.request_ref_scaffolding_files()
        if self.stage_name == "polishing":
            results_list += self.request_polishing_files()
        if self.stage_name == "gap_closing":
            results_list += self.request_gap_closing_files()
        if self.stage_name == "draft_qc":
            results_list += self.request_draft_qc_files()

        return results_list

    def request_common_fasta_stage_files(self):
        # request files common  for fasta generating stages
        results_list = []

        for parameters_label in self.parameters:
            haplotype_list = self.parameters[parameters_label]["haplotype_list"]
            results_list += expand(config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}{extension}",
                                    genome_prefix=[self.config["genome_prefix"],],
                                    assembly_stage=[self.stage_name,],
                                    haplotype=haplotype_list,
                                    extension=[".fasta", ".len", ],
                                    parameters=[parameters_label])

        results_list += self.request_assembly_qc_files()

        return results_list

    def request_draft_qc_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        return results_list

    def request_polishing_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        return results_list

    def request_gap_closing_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        return results_list

    def request_ref_scaffolding_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        return results_list

    def request_dedup_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []
        """
            if config["stage_list"][stage_index] == "purge_dups":
        
        parameters_list = list(stage_dict[current_stage]["parameters"].keys())
        results_list += [
                        *[expand(config["out_dir"] /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label]) for parameters_label in parameters_list],
                        expand(config["out_dir"] / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage],),
                        expand(config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/purge_dups/after.comparison.coverage.png",
                            assembly_stage=[current_stage],
                            parameters=parameters_list
                               ),

                        ]

        for parameters_label in parameters_list:
            if "skipped" not in parameters_label:
                results_list += [[expand(config["out_dir"] / "purge_dups/{parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.{artefact}.fasta",
                                        purge_stage=["first_stage", ] if haplotype == "hap0" else ["first_stage", "second_stage"],
                                        genome_prefix=[config["genome_prefix"], ],
                                        artefact=["junk", "repeat", "haplotig", "ovlp", "highcov"],
                                        haplotype=[haplotype],
                                        parameters=[parameters_label]) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]],
                                 expand(config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
                                     assembly_stage=["purge_dups"],
                                     parameters=[parameters_label]
                                        ),
                                 [expand(config["out_dir"] /  "purge_dups/{parameters}/assembly_qc/purge_dups/{haplotype}/{haplotype}.before-after.comparison.coverage.png",
                                        parameters=[parameters_label],
                                        haplotype=[haplotype],
                                        ) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]]
                                 ]
        """
        return results_list

    def request_hic_scaffolding_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        for parameters_label in self.parameters:
            if "yahs" in parameters_label :
                if not self.config["skip_yahs_hic_file"]:
                    results_list += [expand(config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            assembly_stage=["hic_scaffolding", ],
                                            haplotype=self.parameters[parameters_label]["haplotype_list"],
                                            parameters=[parameters_label])
                                     ]

        return results_list

    def request_hic_alignment_files(self):
        # request specific for this stage files.
        results_list = []
        #for parameters_label in self.parameters:
        #    haplotype_list = self.parameters[parameters_label]["haplotype_list"]
        #    results_list += [expand(config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{haplotype}.{phasing_kmer_length}.rmdup.bam",
        #                            genome_prefix=[self.config["genome_prefix"],],
        #                            assembly_stage=[self.stage_name,],
        #                            haplotype=haplotype_list,
        #                            phasing_kmer_length=[stage_dict[self.stage_name].parameters[parameters_label]["option_set"]["phasing_kmer_length"]],
        #                            parameters=[parameters_label]
        #                            )]
        #
        return results_list

    def request_contig_files(self):
        # request specific for this stage files. fasta, len and qc files are requested by self.request_common_fasta_stage_files method.
        results_list = []

        for parameters_label in self.parameters:
            haplotype_list = self.parameters[parameters_label]["haplotype_list"]

            #TODO: commented code below caused some issues. Check if it is possible to somehow involve alt haplotype or remove this fragment completely
            #if self.parameters[parameters_label]["tool"] == "hifiasm":
            #    if self.parameters[parameters_label]["option_set"]["assembly_ploidy"] > 1:
            #        haplotype_list.append("alt")
            #    else:
            #        haplotype_list.append("alt0")
            #request basic files
            if "hifiasm" in parameters_label:
                results_list += expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}{extension}",
                                        genome_prefix=[self.config["genome_prefix"],],
                                        assembly_stage=["contig",],
                                        haplotype=haplotype_list,
                                        extension=[".unfiltered.gfa.cov", ".unfiltered.gfa.lencov"],
                                        parameters=[parameters_label])

            if self.config["database_set"]["fcs_adaptor"] and (not self.config["skip_fcs_adaptor"]):
                results_list += [expand(config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/contamination_scan/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.report",
                                       genome_prefix=[self.config["genome_prefix"], ],
                                       assembly_stage=[self.stage_name],
                                       haplotype=haplotype_list,
                                       parameters=[parameters_label],
                                       database=self.config["database_set"]["fcs_adaptor"]),
                                ]
            if self.config["database_set"]["fcs"] and (not self.config["skip_fcs"]):
                results_list += [expand(config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/contamination_scan/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.taxonomy",
                                        genome_prefix=[self.config["genome_prefix"], ],
                                        assembly_stage=[self.stage_name],
                                        haplotype=haplotype_list + (["alt" if stage_dict[self.stage_name].parameters[parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if "hifiasm" in parameters_label else []),
                                        parameters=[parameters_label],
                                        database=self.config["database_set"]["fcs"])
                                ]

        results_list += self.request_assembly_qc_files()
        return results_list

    def request_read_phasing_files(self):
        results_list = []
        for datatype in self.config["data_feature_dict"]["phasing"]:
            for parameters_label in list(stage_dict[self.config["phasing_stage"]].parameters.keys()):
                if len(stage_dict[self.config["phasing_stage"]].parameters[parameters_label]["haplotype_list"]) > 1:
                    if datatype in self.config["data_feature_dict"]["paired"]:
                        results_list += [expand(self.config["out_dir"]  / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{pairprefix}{fwd_suffix}{extension}",
                                                fwd_suffix=self.config["data"][datatype]["conv_fwd_sfx"],
                                                extension=self.config["data"][datatype]["conv_ext"],
                                                datatype=[datatype],
                                                stage=[self.config["phasing_stage"], ],
                                                parameters=[parameters_label],
                                                pairprefix=self.config["data"][datatype]["pair_prefix_list"],
                                                genome_prefix=[self.config["genome_prefix"], ],
                                                haplotype=stage_dict[self.config["phasing_stage"]].parameters[parameters_label]["haplotype_list"],
                                                assembly_kmer_length=self.config["assembly_kmer_length"]
                                                ),
                                        ]
                    else:
                        results_list += [expand(self.config["out_dir"]  / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{fileprefix}{extension}",
                                                extension=self.config["data"][datatype]["conv_ext"],
                                                datatype=[datatype],
                                                stage=[self.config["phasing_stage"], ],
                                                parameters=[parameters_label],
                                                fileprefix=self.config["data"][datatype]["conv_file_prefix_list"],
                                                genome_prefix=[self.config["genome_prefix"], ],
                                                haplotype=stage_dict[self.config["phasing_stage"]].parameters[parameters_label]["haplotype_list"],
                                                assembly_kmer_length=self.config["assembly_kmer_length"]
                                                )
                                        ]
        return results_list

    def request_assembly_qc_files(self):
        results_list = []
        if self.config["assembly_qc_level"][self.stage_name] == 0: # skip all the qc
            return results_list

        for parameters_label in self.parameters:
            if self.config["assembly_qc_level"][self.stage_name] >= 1:
                #---- Request stage stat gathering ----
                # it automatically requests Quast and merqury depending on the assembly_qc_level
                # if assembly_qc_level == 1, only quast will be requested,
                # if assembly_qc_level >= 2, both quast and merqury will be requested
                results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                                        genome_prefix=[self.config["genome_prefix"], ],
                                        assembly_stage=[self.stage_name])
                             ]
                #----
                #---- Request microchromosomes ----
                if ("microchromosomes" in self.config) and self.config["microchromosomes"]:
                    results_list += [expand(self.config["out_dir"]  / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.reordered.candidates.microchromosomes.filtered.tsv",
                                            assembly_stage=[self.stage_name],
                                            parameters=[parameters_label],
                                            genome_prefix=[self.config["genome_prefix"],])
                                     ]
                #----

            if self.config["assembly_qc_level"][self.stage_name] >= 3:
                #---- Request assembly based tracks ----
                if not self.config["skip_assembly_tracks"]:
                    assembly_based_track_list = ["gc"]
                    if not self.config["skip_repeats"]:
                        if not self.config["skip_windowmasker"]:
                            assembly_based_track_list.append("windowmasker")
                        if not self.config["skip_trf"]:
                            assembly_based_track_list.append("trf")
                    for track_type in assembly_based_track_list:
                        for window_settings in self.config["qc_settings"]["windows_sets"]:
                            results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                                                   scaffold_length=self.config["qc_settings"]["assembly_scaffold_sets"],
                                                   threshold_type=self.config["qc_settings"]["threshold_types"],
                                                   genome_prefix=[self.config["genome_prefix"], ],
                                                   assembly_stage=[self.stage_name, ],
                                                   track_type=[track_type],
                                                   window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                                                   step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                                                   haplotype=self.parameters[parameters_label]["haplotype_list"],
                                                   parameters=[parameters_label])
                                            ]
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.gap.track.bedgraph",
                                            assembly_stage=[self.stage_name, ],
                                            parameters=[parameters_label],
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            haplotype=self.parameters[parameters_label]["haplotype_list"],)
                                     ]

                if not self.config["skip_telomere"]:
                    if not self.config["skip_telomere_container"]:
                        results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                                        genome_prefix=[self.config["genome_prefix"], ],
                                        assembly_stage=[self.stage_name, ],
                                        parameters=[parameters_label],
                                        haplotype=self.parameters[parameters_label]["haplotype_list"],)
                                       ]
                    if not self.config["skip_telomere_tidk"]:
                        results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_tidk_telomere_all.win{window}.step{window}.track.bedgraph",
                                         window=[parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"]],
                                         genome_prefix=[self.config["genome_prefix"], ],
                                         assembly_stage=[self.stage_name, ],
                                         parameters=[parameters_label],
                                         haplotype=self.parameters[parameters_label]["haplotype_list"],)
                                        ]
                #----
            if self.config["assembly_qc_level"][self.stage_name] >= 4:
                #---- Request WGA ----
                if not self.config["skip_wga"]:
                    if (f"skip_{self.stage_name}_wga" not in self.config) or (not self.config[f"skip_{self.stage_name}_wga"]):
                        results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/wga/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                             query_length=self.config["qc_settings"]["assembly_scaffold_sets"],
                                             target_length=self.config["qc_settings"]["assembly_scaffold_sets"],
                                             genome_prefix=[self.config["genome_prefix"], ],
                                             assembly_stage=[self.stage_name ],
                                             parameters=[parameters_label],
                                             min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                             query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                                 genome_prefix=[self.config["genome_prefix"], ],
                                                                 assembly_stage=[self.stage_name, ],
                                                                 haplotype=self.parameters[parameters_label]["haplotype_list"]),
                                             target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                                 genome_prefix=[self.config["genome_prefix"], ],
                                                                 assembly_stage=[self.stage_name, ],
                                                                 haplotype=self.parameters[parameters_label]["haplotype_list"]),
                                           ),
                                     ]
                        if "reference" in self.config["data"]:
                            results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/wga/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                                     query_length=self.config["qc_settings"]["reference_scaffold_sets"],
                                                     target_length=self.config["qc_settings"]["assembly_scaffold_sets"],
                                                     genome_prefix=[self.config["genome_prefix"], ],
                                                     assembly_stage=[self.stage_name, ],
                                                     parameters=[parameters_label],
                                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                                     query_prefix=list(self.config["data"]["reference"]["ref_dict"].keys()),
                                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                                         genome_prefix=[self.config["genome_prefix"], ],
                                                                         assembly_stage=[self.stage_name, ],
                                                                         haplotype=self.parameters[parameters_label]["haplotype_list"]),)
                                             ]
                #----
                #---- Request ragtag_qc ----
                if not self.config["skip_ragtag_qc"]:
                    if "reference" in self.config["data"]:
                        results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/ragtag/{reference}/{genome_prefix}.{assembly_stage}.{haplotype}.to.{reference}.fasta",
                                                genome_prefix=[self.config["genome_prefix"], ],
                                                assembly_stage=[self.stage_name],
                                                parameters=[parameters_label],
                                                reference=list(self.config["data"]["reference"]["ref_dict"].keys()),
                                                haplotype=self.parameters[parameters_label]["haplotype_list"],
                                                ),
                                         ]
                #----
            if (not self.config["skip_busco"]) and (self.config["assembly_qc_level"][self.stage_name] >= 5):
                #---- Request BUSCO ----
                if not self.config["skip_busco"]:
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}.{busco_version}.tar.gz",
                                            busco_lineage=self.config["busco_lineage_list"],
                                            busco_version=["busco5"],
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            assembly_stage=[self.stage_name, ],
                                            haplotype=self.parameters[parameters_label]["haplotype_list"],
                                            parameters=[parameters_label]),
                                     expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.merged.tsv",
                                            busco_lineage=self.config["busco_lineage_list"],
                                            busco_version=["busco5"],
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            assembly_stage=[self.stage_name],
                                            parameters=[parameters_label]),
                                     expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.merged.tsv",
                                            busco_lineage=self.config["busco_lineage_list"],
                                            busco_version=["busco5"],
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            assembly_stage=[self.stage_name],
                                            haplotype=self.parameters[parameters_label]["haplotype_list"],
                                            parameters=[parameters_label]),
                                     expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.merged.tsv",
                                            busco_lineage=self.config["busco_lineage_list"],
                                            busco_version=["busco5"],
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            assembly_stage=[self.stage_name],
                                            parameters=[parameters_label])
                                     ]
                #----

            if self.config["assembly_qc_level"][self.stage_name] >= 6:
                #---- Request coverage tracks ----
                if not self.config["skip_coverage_tracks"]:
                    for window_settings in self.config["qc_settings"]["windows_sets"]:
                        results_list += [expand(self.config["out_dir"]  / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.coverage_{settings}.{scaffold_length}.win{window}.step{step}.png",
                                                  settings=parameters["tool_options"]["mosdepth"]["options"],
                                                  scaffold_length=self.config["qc_settings"]["assembly_scaffold_sets"],
                                                  window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_settings]["window"],
                                                  step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_settings]["step"],
                                                  genome_prefix=[self.config["genome_prefix"], ],
                                                  assembly_stage=[self.stage_name, ],
                                                  datatype=set(parameters["tool_options"]["assembly_qc"]["coverage"]["datatype_list"]) & set(self.config["data"]),
                                                  haplotype=self.parameters[parameters_label]["haplotype_list"],
                                                  parameters=[parameters_label]),
                                         ]
                #----
                #---- Request purge_dups qc tracks ----
                if not self.config["skip_purge_dups_qc"]:
                    results_list += [expand(self.config["out_dir"]  / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/purge_dups/{datatype}/{genome_prefix}.{assembly_stage}.{haplotype}.dups.extended.bed",
                                            assembly_stage=[self.stage_name,],
                                            parameters=[parameters_label],
                                            haplotype=self.parameters[parameters_label]["haplotype_list"],
                                            datatype=parameters["tool_options"]["assembly_qc"]["purge_dups"]["datatype_list"] ,
                                            genome_prefix=[self.config["genome_prefix"], ],
                                            ) ]
                #----

            if self.config["assembly_qc_level"][self.stage_name] >= 7:
                if not self.config["skip_higlass_mcool"]:
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.{haplotype}.NA.rmdup.higlass.mcool",
                                              haplotype=["reordered" if ("microchromosomes" in self.config) and self.config["microchromosomes"] else "combined"],
                                              genome_prefix=[self.config["genome_prefix"], ],
                                              assembly_stage=[self.stage_name],
                                              parameters=[parameters_label,],)
                                     ]
                if not self.config["skip_hic_for_combined_haplotype"]:
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.{haplotype}.NA.rmdup.pre.mapq{min_mapq}.hic",
                                              min_mapq=[0],
                                              haplotype=["reordered" if ("microchromosomes" in self.config) and self.config["microchromosomes"] else "combined"],
                                              genome_prefix=[self.config["genome_prefix"], ],
                                              assembly_stage=[self.stage_name],
                                              parameters=[parameters_label,],)
                                     ]

                if not self.config["skip_pretext"]:
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.pretext.track.info",
                                              haplotype=["reordered" if ("microchromosomes" in self.config) and self.config["microchromosomes"] else "combined"],
                                              genome_prefix=[self.config["genome_prefix"], ],
                                              assembly_stage=[self.stage_name],
                                              parameters=[parameters_label,],)
                                     ]
                    # request pretext map for a whole genome
                    results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.{haplotype}.NA.{subset}.rmdup.mapq{mapq}.{res}.tracks.pretext",
                                              res=parameters["tool_options"]["pretextmap"]["res"],
                                              haplotype=["reordered" if ("microchromosomes" in self.config) and self.config["microchromosomes"] else "combined"],
                                              subset=["all"],
                                              genome_prefix=[self.config["genome_prefix"], ],
                                              assembly_stage=[self.stage_name],
                                              parameters=[parameters_label,],
                                              mapq=parameters["tool_options"]["pretextmap"]["mapq"],)
                                     ]

                    if candidate_chr_id_list:
                        # request pretext map for curation units (if provided)
                        results_list += [expand(self.config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/NA/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.NA.{subset}.rmdup.precurated.mapq{mapq}.{res}.tracks.pretext",
                                              res=parameters["tool_options"]["pretextmap"]["res"],
                                              haplotype=["reordered" if ("microchromosomes" in self.config) and self.config["microchromosomes"] else "combined"],
                                              subset=candidate_chr_id_list,
                                              genome_prefix=[self.config["genome_prefix"], ],
                                              assembly_stage=[self.stage_name],
                                              parameters=[parameters_label,],
                                              mapq=parameters["tool_options"]["pretextmap"]["mapq"])
                                     ]
        return results_list

    def request_filter_reads_files(self):
        results_list = []

        for datatype in self.config["data_feature_dict"]["filter"]:
            results_list += [expand(self.config["out_dir"] / "data/{datatype}/filtered/{fileprefix}{extension}",
                                    datatype=[datatype,],
                                    extension=[self.config["data"][datatype]["conv_ext"]],
                                    fileprefix=self.config["data"][datatype]["conv_file_prefix_list"])]

        return results_list

    def request_mtdna_files(self):
        # request specific for this stage files.
        results_list = []
        if not self.config["skip_mitohifi_reads"]:
            if not self.config["skip_mitohifi_reads_per_file"]:
                if "hifi" in self.config["data"]:
                    results_list += [ expand(self.config["out_dir"] / "mtDNA/mitohifi/{mtdna_ref}/hifi/final/{fileprefix}/FINISH_FLAG",
                                             mtdna_ref=["recommended"],
                                             fileprefix=self.config["data"]["hifi"]["conv_file_prefix_list"])]
            if not self.config["skip_mitohifi_reads_combined"]:
                if "hifi" in self.config["data"]:
                    results_list += [expand(self.config["out_dir"] / "mtDNA/mitohifi/{mtdna_ref}/hifi/combined/hifi.combined/FINISH_FLAG",
                                            mtdna_ref=["recommended"],)]
        if not self.config["skip_mitoz"]:
            if (not self.config["skip_mitoz_hic"]) and ("hic" in self.config["data"]):
                results_list += [expand(self.config["out_dir"] / "mtDNA/mitoz/denovo/{datatype}/{stage}/{pairprefix}/FINISH_FLAG",
                                        datatype=["hic",],
                                        stage=["final"],
                                        pairprefix=self.config["data"]["hic"]["pair_prefix_list"])]
            if (not self.config["skip_mitoz_illumina"]) and ("illumina" in self.config["data"]):
                results_list += [expand(self.config["out_dir"] / "mtDNA/mitoz/denovo/{datatype}/{stage}/{pairprefix}/FINISH_FLAG",
                                        datatype=["illumina",],
                                        stage=["final"],
                                        pairprefix=self.config["data"]["illumina"]["pair_prefix_list"])]

        return results_list

    def request_kmer_qc_files(self, stage):
        results_list = []
        for datatype in self.config["data_feature_dict"]["genome_size"]:
            if (datatype == "hic") and (self.config["skip_hic_genomescope"]):
                continue
            for kmer_tool in parameters["tool_options"]["kmer_qc"]["kmer_counter_list"]:
                results_list += [expand(self.config["out_dir"] / "kmer/{datatype}/{stage}/{analysis_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{analysis_tool}.parameters",
                                       datatype=[datatype,],
                                       genome_prefix=[self.config["genome_prefix"], ],
                                       analysis_tool=parameters["tool_options"]["kmer_qc"]["genome_size_estimators"],
                                       stage=[stage,],
                                       kmer_tool=[kmer_tool,],
                                       kmer_length=parameters["tool_options"][kmer_tool][datatype]["kmer_length"],
                                     )]

            if not self.config["skip_per_lib_genome_estimation"]: # per lib estimation is possible only for meryl kmer counter, as other for other kmer counter per-lib databases are not calculated
                kmer_tool = "meryl"
                results_list += [expand(self.config["out_dir"] / "kmer/{datatype}/{stage}/{analysis_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.{analysis_tool}.parameters",
                                           datatype=[datatype,],
                                           genome_prefix=[self.config["genome_prefix"], ],
                                           analysis_tool=parameters["tool_options"]["kmer_qc"]["genome_size_estimators"],
                                           stage=[stage,],
                                           kmer_tool=[kmer_tool,],
                                           read_prefix=self.config["data"][datatype]["pair_prefix_list"] if datatype in self.config["data_feature_dict"]["paired"] else self.config["data"][datatype]["conv_file_prefix_list"],
                                           kmer_length=parameters["tool_options"][kmer_tool][datatype]["kmer_length"],
                                        )]
        # TODO: issues with draw_gc_plot.py script from KrATER, fix it later
        """ 
        if not self.config["skip_kmer_gcp"]:
            for datatype in set(parameters["tool_options"]["gcp"]) & set(self.config["data"]):
                results_list += [expand(self.config["out_dir"]/ "kmer/{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.heatmap.png",
                                        datatype=[datatype,],
                                        stage=["final",],
                                        kmer_length=parameters["tool_options"]["gcp"][datatype]["kmer_length"],
                                        min_coverage=parameters["tool_options"]["gcp"][datatype]["min_coverage"],
                                       )]
        """


        if not self.config["skip_kmer_smudgeplot"]:
            for datatype in self.config["data_feature_dict"]["genome_size"]:
                for kmer_tool in parameters["tool_options"]["kmer_qc"]["kmer_counter_list"]:
                    results_list += [expand(self.config["out_dir"]/ "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_warnings.txt",
                                               lower_boundary=parameters["tool_options"]["smudgeplot"][datatype]["lower_boundary"],
                                               upper_boundary=parameters["tool_options"]["smudgeplot"][datatype]["upper_boundary"],
                                               datatype=[datatype,],
                                               stage=["final",],
                                               kmer_tool=[kmer_tool,],  
                                               kmer_length=parameters["tool_options"][kmer_tool][datatype]["kmer_length"],
                                               ),
                                    expand(config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.smudgeplot.boundaries",
                                              datatype=[datatype,],
                                              stage=["final",],
                                              kmer_tool=[kmer_tool,],
                                              kmer_length=parameters["tool_options"][kmer_tool][datatype]["kmer_length"],
                                              )
                                    ]

        #----



        return results_list

    def request_read_contamination_scan_files(self):
        results_list = []
        for datatype in self.config["data_feature_dict"]["kraken"]:
            results_list += [expand(config["out_dir"] / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
                                   datatype=[datatype],
                                   database=config["database_set"]["kraken2"],
                                   )
                            ]
        return results_list

    def request_qc_files(self, stage):
        results_list = []
        if not self.config["skip_fastqc"]:
            results_list += [[expand(self.config["out_dir"] / "qc/fastqc/{fastqc_datatype}/{stage}/{fileprefix}_fastqc.zip",
                               fastqc_datatype=[dat_type, ],
                               stage=[stage, ],
                               fileprefix=self.config["data"][dat_type]["conv_file_prefix_list"]) for dat_type in self.config["data_feature_dict"]["fastqc"]],
                             expand(self.config["out_dir"] / "qc/multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=self.config["data_feature_dict"]["fastqc"],
                             stage=[stage,]),]
        if not self.config["skip_nanoplot"]:
            results_list += [expand(self.config["out_dir"] / "qc/nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                               datatype=self.config["data_feature_dict"]["long_read"],
                               stage=[stage, ],
                               )]
        if not self.config["skip_nanoqc"]:
             results_list += [[expand(self.config["out_dir"] / "qc/nanoqc/{datatype}/{stage}/{fileprefix}",
                               datatype=[dat_type, ],
                               stage=[stage, ],
                               fileprefix=self.config["data"][dat_type]["conv_file_prefix_list"]) for dat_type in self.config["data_feature_dict"]["long_read"]]
                             ]
        if not self.config["skip_tadbit"]:
            if ("hic" in self.config["data"]) and ((self.config["hic_enzyme_set"] == "custom") or self.config["hic_enzyme_dict"][self.config["hic_enzyme_set"]]):
                results_list += [expand(self.config["out_dir"] / "qc/tadbit/hic/{stage}/{genome_prefix}.tadbit.stats",
                                 stage=[stage, ],
                                 genome_prefix=[self.config["genome_prefix"]])]

        return results_list

    def parse_stage_parameters_from_path(self, path):
        string_list = str(path).split("/") # in case if string is nota str, but a Path

        candidate_stage_list = []
        candidate_stage_index_list = []
        for allowed_stage in self.config["allowed_stage_list"]:
            for index in range(0, len(string_list)):
                if string_list[index] == allowed_stage:
                    candidate_stage_list.append(allowed_stage)
                    candidate_stage_index_list.append(index)
        if len(candidate_stage_list) > 1:
            raise ValueError(f"ERROR!!! Impossible to detect a stage for path: {path}. Candidate stages: {', '.join(candidate_stage_list)} ")
        elif len(candidate_stage_list) == 0:
            raise ValueError(f"ERROR!!! Impossible to detect a stage for path: {path}. No allowed stage was found in it.")
        else:
            stage = candidate_stage_list[0]
            element_index = candidate_stage_index_list[0]

        parameters = string_list[element_index + 1]

        prev_parameters = None if ".." not in parameters else "..".join(parameters.split("..")[:-1])

        result_dict = {"stage": stage,
                        "parameters": parameters,
                        "prev_parameters": prev_parameters}

        return result_dict

    @staticmethod
    def get_prev_stage_parameters(parameters):
        return "..".join(parameters.split("..")[:-1])
