tmux attach -t atlas_transfer 2>/dev/null || tmux new -s atlas_transfer

project="LCSB_CBG"

smb-storage mount guillem.santamaria --project "$project"


project_lower="${project,,}"

proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev"

dir_base=$(basename $proj_dir)

target_dir=~/atlas.uni.lux-${project_lower}/Guillem/

mkdir -p $target_dir

rsync -av --progress --whole-file --inplace --dry-run \
        --no-perms --no-owner --no-group \
        $proj_dir "$target_dir"

rsync -av --progress --whole-file --inplace \
        --no-perms --no-owner --no-group \
        $proj_dir "$target_dir"