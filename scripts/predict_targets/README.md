# Intructions for launching training scripts

## Training scripts with shuffle controls
This script launches simultaneosly a training script working with the real labels, and a job with shuffled labels. Then, it calls compare_runs.py to obtain a comparison of the results. It has different working modes.


```bash
bash submit_with_shuffle_control.sh tf                   # Submits real + shuffled + comparison (on tf data)
bash submit_with_shuffle_control.sh expr --real-job <id> # Submits real shuffle + comparison (on expression data). In case the real job is already running.
bash submit_with_shuffle_control.sh mc   --real-job <id> # Submits real shuffle + comparison (on MetChange data). In case the real job is already running.
```

