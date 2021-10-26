#!/bin/bash


mkdir res
mkdir res/uc
for graph in  dtelekom; do
  for cache in TBGRD; do
    if [ $cache == 'TBGRD' ]; then
      for freq in -1; do
        for param in .0005 ; do
          for seed in {51..81}; do
            echo python CacheNetwork.py res/uc/uc \
              --graph_type $graph \
              --cache_type $cache \
              --graph_size 100 \
              --min_capacity 1 \
              --max_capacity 5 \
              --query_nodes 5 \
              --catalog_size 50 \
              --demand_size 200 \
              --max_weight 1000 \
              --time 1000 \
              --action_selector_eta $param \
              --correlated_action_selectors True \
              --random_seed $seed \
              --frequency $freq \
              --trace_location traces/changing_popularity_catalog_50.pkl &
          done
        done
      done
    fi
  done
done
sleep 1000
