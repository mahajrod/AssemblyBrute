localrules: create_local_links

rule create_local_links: # abstract rule to use via redefining
    input:
        object="filename",
        log_dir=ancient("{link_dir}/log")
    output:
        object="{link_dir}/{link_name}filename"
    log:
        ln="{link_dir}/log/create_local_link.{link_name}.ln.log",
    resources:
        queue=config["queue"]["cpu"]["name"],
        cpus=parameters["threads"]["create_links"],
        node_options=parse_node_list("create_links"),
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"],
    threads:
        parameters["threads"]["create_links"]
    run:
        input_dict = dict(input.items())
        input_dict.pop("log_dir")
        output_dict = dict(output.items())

        with open(log.ln, "w") as log_fd:
            with redirect_stdout(log_fd), redirect_stderr(log_fd):
                for entry in input_dict:
                    create_relative_link(input_dict[entry], output_dict[entry])
