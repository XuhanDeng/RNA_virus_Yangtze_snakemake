# AlphaFold3 install + database download workflow, Singularity-based.
# https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md
#
# DelftBlue has no Docker daemon (see
# https://doc.dhpc.tudelft.nl/delftblue/howtos/singularity/), so the
# "official" build path (docker build -> local registry -> singularity build
# docker://localhost...) doesn't work on the cluster itself.
#
# `apptainer build --fakeroot` directly from AF3's Dockerfile was tried and
# does NOT work here either: fakeroot itself is fine on this account, but
# Apptainer's builder does not understand plain Dockerfile syntax (it only
# accepts its own .def format, or docker://.../oras://... image references),
# so it fails immediately with "invalid header keyword found: from ...".
#
# So this pulls the prebuilt community image instead:
#   https://github.com/colbyford/af3-on-hpc
#   https://hub.docker.com/r/cford38/alphafold3
# NB: that image is pinned at AF3 v3.0.0 (never updated past that per the
# af3-on-hpc README, which itself recommends staying on v3.0.0 due to import
# errors reported in v3.0.1). Use this to get running now; if a current
# version is needed later, the only paths are (a) build on a machine that has
# Docker and scp the .sif to config["alphafold3"]["sif"], or (b) hand-write
# an Apptainer .def file translating the Dockerfile's stages.
#
# Model weights (af3.bin.zst) are documented as gated, but
# config["alphafold3"]["weights_url"] currently serves the file with no auth
# challenge (verified 2026-09-01: HTTP 200, ~973 MiB) -- downloaded
# automatically below. Usage of the weights is still governed by
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md
# regardless of download access; if that URL ever starts requiring auth,
# download manually and place af3.bin.zst under config["alphafold3"]["model_dir"].

configfile: "config/config.yaml"

_AF3 = config["alphafold3"]


rule all:
    input:
        _AF3["sif"],
        _AF3["db_marker"],
        _AF3["model_dir"] + "/af3.bin.zst",


rule clone_alphafold3:
    output:
        marker = _AF3["install_dir"] + "/.clone.done",
    log:
        log = "log/00_alphafold3/clone_alphafold3.log",
        err = "log/00_alphafold3/clone_alphafold3.err",
    params:
        install_dir = _AF3["install_dir"],
        repo_url    = _AF3["repo_url"],
        version     = _AF3["version"],
    shell:
        """
        mkdir -p $(dirname {log.log})

        if [ ! -f {params.install_dir}/docker/Dockerfile ]; then
            rm -rf {params.install_dir}
            git clone --branch {params.version} --depth 1 \
                {params.repo_url} {params.install_dir} \
                > {log.log} 2> {log.err}
        else
            echo "Repo already cloned, skipping." > {log.log}
        fi

        touch {output.marker}
        """


# Pull the prebuilt community image directly -- no Docker daemon, no
# fakeroot/Dockerfile-parsing needed. Pinned at AF3 v3.0.0 (see module
# docstring above for why); swap config["alphafold3"]["docker_image"] or
# replace this rule's output file if/when a newer prebuilt image or a
# Docker-built .sif becomes available.
rule build_alphafold3_sif:
    output:
        sif = _AF3["sif"],
    log:
        log = "log/00_alphafold3/build_sif.log",
        err = "log/00_alphafold3/build_sif.err",
    params:
        image = _AF3["docker_image"],
    shell:
        """
        mkdir -p $(dirname {log.log}) $(dirname {output.sif})
        singularity pull {output.sif} docker://{params.image} \
            > {log.log} 2> {log.err}
        """


# Genetic databases (BFD-small, MGnify, PDB mmCIF, PDB seqres, UniProt,
# UniRef90, NT, Rfam, RNACentral): ~252 GB download, ~630 GB uncompressed.
rule download_alphafold3_dbs:
    input:
        marker = _AF3["install_dir"] + "/.clone.done",
    output:
        marker = _AF3["db_marker"],
    log:
        log = "log/00_alphafold3/download_dbs.log",
        err = "log/00_alphafold3/download_dbs.err",
    conda:
        "../envs/alphafold3.yaml"
    params:
        install_dir = _AF3["install_dir"],
        db_dir      = _AF3["db_dir"],
    shell:
        """
        mkdir -p $(dirname {log.log}) {params.db_dir}
        db_abs=$(realpath {params.db_dir})

        python {params.install_dir}/fetch_databases.py \
            --download_destination $db_abs \
            > $(realpath {log.log}) 2> $(realpath {log.err})

        touch {output.marker}
        """


# Model parameters (~973 MiB). See module docstring re: gated-but-currently-
# public access; usage remains subject to WEIGHTS_TERMS_OF_USE.md.
rule download_alphafold3_weights:
    output:
        weights = _AF3["model_dir"] + "/af3.bin.zst",
    log:
        log = "log/00_alphafold3/download_weights.log",
        err = "log/00_alphafold3/download_weights.err",
    params:
        url = _AF3["weights_url"],
    shell:
        """
        mkdir -p $(dirname {log.log}) $(dirname {output.weights})
        wget -c {params.url} -O {output.weights} \
            > {log.log} 2> {log.err}
        """
