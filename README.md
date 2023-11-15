# RNA-Seq Pipeline

## Downloading and Installing the code

1. Setting up ssh on the Broad server

    Follow the instructions below to create an SSH key and connecting it to your GitHub account

    https://docs.github.com/en/authentication/connecting-to-github-with-ssh

2. Download the code on to your $HOME folder in the Broad server
    ```bash
    # go to your $HOME folder
    cd 
    # download the code
    git clone git@github.com:hsethura/RNA-Seq-Pipeline.git
    ```

## Guidelines
### Acquiring a node on the Broad server
```bash
use UGER
ish
```

### Configuration files
There are two configuration files that are currently in use:

* ```idi/moc_ec/MOC/config_files/PC_config.yaml``` for RtS projects
* ```idi/moc_ec/MOC/config_files/PC_config_BacDrop.yaml``` for SCR projects

### Running commands
It is advised to run all commands from the project root directory, i.e., 

```bash
cd $HOME/RNA-Seq-Pipeline
```

### Relative paths

All relative paths mentioned below are with respect to $HOME/RNA-Seq-Pipeline

## Pipeline commands

### Add URL to the “WB_linking” sheet

Following completion of the MOCS-XXXX_Pool_Sub_WB, add it’s URL to the “WB_linking” sheet of [RtS MOCS DB](https://docs.google.com/spreadsheets/d/1WC8s_Y6uMbOCXNCQMIBWxBxglHOIo0HHlH5G51U7e1U/edit#gid=426292915) and insure the pools appear in the “MOCS DB” sheet.

### Moving data from walkup moc server

Once email is received indicating sequencing has been completed for a flowcell, move data from walkup moc server.  For this you will need the MOCS-ID and the flowcell ID, the latter can be found in the walkup email subject line.  Run the following command

```bash
sh idi/moc_ec/MOC/scripts/MOC_wu_move.sh <MOCS-ID> <FLOWCELL-ID>
```

**Example:**
```
sh /idi/moc_ec/MOC/scripts/MOC_wu_move.sh MOCS-0004 H52G2DMXY	
```

**Command line options:**
```
-move_data: include if you want to move data from getsite.  (Default Y)
-conf: path to config file.  (Default /idi/moc_ec/MOC/config_files/Universal_config.yaml)
-metrics: run metrics script (Default Y)
-raw_seq_path: the target location for the data (Defaults path designated by “Seq_base:” in the config file)
```

All files will be generated in a directory named MOCSID (with hyphen removed) in the path designated by “Seq_base:” in the config file (/idi/moc_ec/RawSeq_data/RtS/ in default config) or -raw_seq_path 

Metrics file will be emailed to you and copied into the seq dir.

> **Note:** Usually there must be the same number of files per lane in the metrics file. If there are uneven read counts among each lane it may mean an incomplete download.  Repeat and review the submission (i.e. check that there weren’t different pools submitted for each lane). Contact Tammy (tamason@broadinstitute.org) and report the error.

**Example of an incomplete download:**
```FC: H2527BGX9
LANE 1: 168
LANE 2: 168
LANE 3: 78
LANE 4: 42
```

**Example of a complete download:**
```
FC: H2J32BGX9
LANE 1: 168
LANE 2: 168
LANE 3: 168
LANE 4: 168
```


This script runs another script that sets symlinks for all of the projects.  If there are changes to the MOCS spreadsheet this needs to be rerun to correct the symlinks

```
sh idi/moc_ec/MOC/scripts/MOCS_move.sh -limit <MOCS_ID>
```

### Review the metrics file. 

Verify that all flowcells and all lanes generated sufficient read numbers:

**Expected output:**
```
NextSeq: >400M reads total, >80%matched to an index, even coverage across 4 lanes
```

If total reads or “orphan” rates are significantly lower or higher, respectively than indicated above, look to see if there are a high % of unmatched reads aligning to known index sequences.  If the data are dual barcoded, run the following script to look for P7 and P5 pairs enriched in the unmatched reads:

```
sh idi/moc_ec/MOC/scripts/dual_orphan_find.sh <path to *1.unmatched.barcode_1.fastq.gz>
```

Contact Jonathan to get input on possible reasons for these discrepancies and discuss the possible need to contact the walkup team for resequencing/redemultiplexing.

### Copy paste file contents into “MOCS metrics”

Copy paste file contents into “MOCS metrics” sheet of MOCS WB (you may need to “split text to column” in Data scroll down menu) to compare output to target depth

### Launching analysis pipeline

Launch analysis pipeline with the following command:
```
sh idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh <MOC-ID> -user_id [script options] :[pipeline options]:
```

**Example:**
```
sh idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh MOCP-0035 -user_id 

sh idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh MOCP-0001 -user_id -move_key N :--no_split --no_merge --no_align:
```

**Running the pipeline with host analysis:**
```
sh idi/moc_ec/MOC/scripts/MOC_RtS_pipe_v2.sh DRS_0001.3 -user_id :--do-host:
```

**Script options (see script for full list)**
```
# -conf: sets path to config file (default idi/moc_ec/MOC/config_files/Universal_config.yaml)
# -q: sets name of header for Q_VAL (default MOC-ID)
# -move_key: setting to N skips moving google sheet to server (default Y)
# -use_p7: indicates fastq file is demultiplexed by P7 index (default Y)
# -use_p5: indicates fastq file is demulitplexed by P5 index (default N)
# -no_pipe: including or setting to Y skips running pipeline (default N)
# -qsub: run pipeline in UGER (default N)
# -gid: enter google ID if missing from DB (default 0 => get it from DB)
# -smocid: enter smocid if missing from DB (default 0 => get it from DB)
# -import_gs: updates DBs from google drive on server (default Y)
# -user_id: add userID to temp and results paths (default N.  If Y, add login name, if value added after option that userID will be added)
# -moc_id: add MOC-ID to temp and results paths (default Y)
# -move_ref: move reference files from Gdrive/server to ref directory  (default Y)
# -raw_seq_path: path to raw data (overrides path in config file) * remember to include the terminal “/”
# -symlink: path to directory where symlinks to raw data are located
# -key_sheet: name of sheet in key file to be used by pipeline (default "Sample_Information")
# -split:  split fastqs by inline barcode (default Y)
# -merge:  merge split fastqs (default Y)
# -align:  align merged fastqs to reference (default Y)
# -count:  count reads per features from aligned bams (default Y)
# -fastqc: fastqc on the raw files. Fastqc files are written to $MOC_SYM_DIR\analysis\fastqc (default Y)
# -fastqc_sample: fastqc on the sample-level files. Fastqc files are written to $MERGE_DIR\analysis\fastqc (default Y)
# -adapter_analysis: Adapter analysis on read 2 of the raw files. Adapter analysis files are written to $MOC_SYM_DIR\analysis\adapter_analysis (default Y)
# -adapter_analysis_sample: Adapter analysis on read 2 of the sample-level files. Adapter analysis files are written to $MERGE_DIR\analysis\adapter_analysis (default Y)
# -trimmomatic: Adapter trimming using Trimmomatic. Trimmomatic files are written to $MOC_SYM_DIR\analysis\trimmomatic (default N)
# -trimmomatic_sample: Adapter trimming using Trimmomatic. Trimmomatic files are written to $MERGE_DIR\analysis\trimmomatic_sample (default Y)
```

**Pipeline options**
The default pipeline options are the following:
```
 --ADD3 30 --ADD5 20 --remove_splitted --do_patho --MOC_id_ref $MOC_ID
 ```

> **Note**: If you’d like to add these options to the default options them enter :<other options> on command line.  If you want to override them with other option enter :<other options --no_def_opt>: MAKE SURE TO INCLUDE “:” ON EITHER SIDE OF THESE OPTIONS. 

One key option is --MOC_id_ref.  Use this followed by the MOC-ID when moving ref files from GDrive using GDrive_fetch.sh.
 
Other commonly used pipeline options are described [here](https://docs.google.com/document/d/1BiWyl7-nz7QPbAwAXIZiKD6HeLVVTNlieBTQCU_nynE/edit).

If you have an error indicating authentication issues (if the word json appears in the error output) follow these [instructions](https://docs.google.com/document/d/16sVTq8PN4Er8E2GUUp2jTWSsalV8IWlhyD0lXtWrJXw/edit)

**Launching analysis pipeline using UGER**

If you want to launch the pipeline through UGER, use the following command:
```
 sh idi/moc_ec/MOC/scripts/QSUB_RtS_analysis_launch.sh <MOC-ID> [script options] :[pipeline options]:
```

Uses same options as regular script.  Creates a directory $RESULTS_PATH”/UGER/” into which the UGER err and out files appended with a timestamp are written.


If you have resequenced some of the samples and want to merge the old and new data, add the information for the redo samples to the key and give the samples their original sample ID appended by “_add”.  Then after downloading the data from both seq runs run this in lieu of the normal pipeline:

```
sh idi/moc_ec/MOC/scripts/redo_merge_RtS.sh <MOC-ID> [script options]
```

> **NOTE**: The sequences of inline barcodes are here:

```
/broad/idi/moc_ec/MOC/files/rts_bcs.txt

/broad/idi/moc_ec/MOC/files/scr_bcs.txt
```
If changed ensure no blank lines at the end of the file

### Running the pipeline to get only split metrics

If you’re just splitting and merging and want split metrics, run the following command:
```
sh idi/moc_ec/MOC/scripts/Split_distribution_metrics.sh MOC-XXX
```

If someone else ran the initial split and merged, add the following option:
-uid <userid of person who ran the pipeline>

### Running DESeq2

As long as CG_IDs and CG_pairs are filled in in the key file, DESeq2 and edgeR will be automatically launched.  If you want to launch these independently use the following command:

```
sh idi/moc_ec/MOC/scripts/Key_DE_launch2.sh <MOC-ID>
```

By default, this will look for files in the results and temp paths in the config file WITHOUT the userID and WITH the MOC ID.  To change this use the -moc_id and -user_id options.  For example:

```
sh idi/moc_ec/MOC/scripts/Key_DE_launch2.sh <MOC-ID> -user_id -moc_id N
```

...will append your userID and not a MOC-ID

## Data Interpretation

A guide for interpreting data from the pipeline can be found [here](https://docs.google.com/document/d/1FXtmN9AWRTp3VoGEgodAb9QddS0Uyz6R/edit)

## Data Transfer

To share data with collaborators, run the following command 

```
sh idi/moc_ec/MOC/scripts/MOC_data_transfer2.sh <MOC-ID>
```
	
**Options:**
```
# -fastq: include fastq files along with bams in Aspera (default N)
# -internal: is this an internal or external project (default whatever is in DB)
# -import_gs: import databases (default N)
# -import_key: import key file (default N)
# -user_id: add if analysis was run with this option.  
```

This should create a directory in Transfer_path (/broad/hptmp/MOC/transfer/MOC) to which the user guide will be added.  All result files, DE directory, and user guide will be moved to a results directory in the MOC-ID directory of the gdrive that contains the key and ref files.  If an internal project, you should get an email you can forward to the collaborator telling them where the data live on the server.  If an external project, you should get an email with an attached manifest and instructions that you can use to set up an Aspera site [here](https://transfer.broadinstitute.org/auth/login).  Once the site is set up you can forward that email to the collaborator. 

Once you’ve sent the email(s), enter the date in the “Data DB” sheet in the [MOC Production DB](https://docs.google.com/spreadsheets/d/1tF0Cc6CwbT_bS3JLKIFGA3oPK9IBqjilyu3gtgmYDas/edit#gid=290518403)

