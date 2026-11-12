#!/bin/bash

# Two versions for matching negatives and forcing the model to learn from actual signatures and not the priors (fixed protein info in the model)

# Submit the version A: plain, add anegatives per protein.

bash submit_with_shuffle_control.sh mc_matched
