script_path="$(readlink -f "$0")"
echo $script_path
config_path="$(dirname $(dirname $script_path))"/config_files/PC_config.yaml
echo $config_path
