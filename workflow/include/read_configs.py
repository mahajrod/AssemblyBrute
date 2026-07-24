

#---- Read config files ----
#-------- Read core config file --------
with open(config["main_config_file"], "r") as core_yaml_fd:
    config.update(yaml.safe_load(core_yaml_fd))
#--------
#-------- Read secondary tools config file -------
with open(config["secondary_tool_config_file"], "r") as secondary_tool_fd:
    copy_absent_entries(yaml.safe_load(secondary_tool_fd), config)
#--------
#-------- Read database and container config file --------
with open(config["database_and_container_config_file"], "r") as databases_and_containers_fd:
    copy_absent_entries(yaml.safe_load(databases_and_containers_fd), config)
#--------
#------- Read configs of coretools with separated configs -------
coretool_config_dir_path = Path(config["coretool_config_dir"])
for parameter_set in config["parameters"]:
    for coretool in config["coretool_config_dict"]:
        if (coretool_config_dir_path / parameter_set).exists():
            config["parameters"][parameter_set]["tool_options"][coretool] = {}
            if (coretool_config_dir_path / parameter_set / f"{coretool}.tab").exists():
                coretool_config_df = pd.read_csv(coretool_config_dir_path / parameter_set / f"{coretool}.tab", sep="\t",
                                                 header=0, index_col=0,
                                                 converters=parse_coretool_config_converters(config["coretool_config_dict"][coretool]),
                                                 dtype=parse_coretool_config_datatypes(config["coretool_config_dict"][coretool]))
                copy_absent_entries(coretool_config_df.to_dict(orient='index'),
                                    config["parameters"][parameter_set]["tool_options"][coretool])
                #print(config["parameters"][parameter_set]["tool_options"][coretool])
#--------

#-------- Read cluster config file ---------
with open(config["cluster_config_file"], "r") as cluster_fd:
    copy_absent_entries(yaml.safe_load(cluster_fd), config)
#--------

#-------- Read datatype config file ---------
with open(config["datatype_config_file"], "r") as cluster_fd:
    copy_absent_entries(yaml.safe_load(cluster_fd), config)
#--------

#-------- Read 'skip' config file --------
with open(config["skip_config_file"], "r") as skip_yaml_fd:
    for key, value in yaml.safe_load(skip_yaml_fd).items():
        if key not in config:
            config[key] = value
#--------
#-------- Read resources config files --------
resources_dir_path = Path(config["resources_dir"])
for resource, res_datatype in zip(["threads", "memory_mb", "time"], [int, int, str]):
    resource_df = pd.read_csv(resources_dir_path / f"{config['resource_profile']}/{resource}.tab",
                              sep="\t", header=0, index_col=0)
    for config_label in resource_df.columns:
        config["parameters"][config_label][resource] = resource_df[config_label].to_dict(into=OrderedDict)

#--------