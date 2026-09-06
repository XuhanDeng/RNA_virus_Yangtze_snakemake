# LucaProt install + database download workflow
# Run on login node — requires internet access.

configfile: "config/config.yaml"


rule all:
    input:
        config["lucaprot"]["marker_install"],
        config["lucaprot"]["marker_db"]


rule install_lucaprot:
    output:
        marker = config["lucaprot"]["marker_install"]
    log:
        log = "log/00_lucaprot/install.log",
        err = "log/00_lucaprot/install.err"
    conda:
        "../envs/lucaprot.yaml"
    params:
        install_dir = config["lucaprot"]["install_dir"],
        repo_url    = config["lucaprot"]["repo_url"]
    shell:
        """
        mkdir -p $(dirname {log.log})

        if [ ! -f {params.install_dir}/src/predict_one_sample.py ]; then
            rm -rf {params.install_dir}
            git clone {params.repo_url} {params.install_dir} \
                > {log.log} 2> {log.err}
        else
            echo "Repo already cloned, skipping." > {log.log}
        fi

        touch {output.marker}
        """


rule download_lucaprot_db:
    input:
        marker_install = config["lucaprot"]["marker_install"]
    output:
        marker = config["lucaprot"]["marker_db"]
    log:
        log = "log/00_lucaprot/download_db.log",
        err = "log/00_lucaprot/download_db.err"
    conda:
        "../envs/lucaprot.yaml"
    params:
        install_dir = config["lucaprot"]["install_dir"],
        db_dir      = config["lucaprot"]["db_dir"],
        emb_dir     = config["lucaprot"]["db_dir"] + "/emb"
    shell:
        """
        mkdir -p $(dirname {log.log})
        mkdir -p {params.emb_dir}
        mkdir -p {params.db_dir}
        log_abs=$(realpath {log.log})
        err_abs=$(realpath {log.err})
        emb_abs=$(realpath {params.emb_dir})
        db_abs=$(realpath {params.db_dir})
        marker_abs=$(realpath {output.marker})

        cd {params.install_dir}
        python src/predict_one_sample.py \
            --protein_id protein_1 \
            --sequence MTTSTAFTGKTLMITGGTGSFGNTVLKHFVHTDLAEIRIFSRDEKKQDDMRHRLQEKSPELADKVRFFIGDVRNLQSVRDAMHGVDYIFHAAALKQVPSCEFFPMEAVRTNVLGTDNVLHAAIDEGVDRVVCLSTDKAAYPINAMGKSKAMMESIIYANARNGAGRTTICCTRYGNVMCSRGSVIPLFIDRIRKGEPLTVTDPNMTRFLMNLDEAVDLVQFAFEHANPGDLFIQKAPASTIGDLAEAVQEVFGRVGTQVIGTRHGEKLYETLMTCEERLRAEDMGDYFRVACDSRDLNYDKFVVNGEVTTMADEAYTSHNTSRLDVAGTVEKIKTAEYVQLALEGREYEAVQ \
            --emb_dir $emb_abs \
            --truncation_seq_length 4096 \
            --dataset_name rdrp_40_extend \
            --dataset_type protein \
            --task_type binary_class \
            --model_type sefn \
            --time_str 20230201140320 \
            --step 100000 \
            --threshold 0.5 \
            --gpu_id -1 \
            --torch_hub_dir $db_abs \
            > $log_abs 2> $err_abs

        touch $marker_abs
        """
