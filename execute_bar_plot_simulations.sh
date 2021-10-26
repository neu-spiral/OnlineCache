#!/bin/bash
#source ~/.bashrc
#conda activate phdwork

path=psn
trace=fixed_popularity_catalog_50.pkl

mkdir res
mkdir res/$path
for graph in balanced_tree erdos_renyi balanced_tree hypercube geant abilene dtelekom; do
  for cache in TBGRD ; do
    if [ $cache == 'EWMAGRAD' ]; then
      for param in .005 .01 .05 .1; do
        echo python CacheNetwork.py res/$path/$path \
          --graph_type $graph \
          --cache_type $cache \
          --graph_size 100 \
          --min_capacity 1 \
          --max_capacity 5 \
          --query_nodes 5 \
          --catalog_size 50 \
          --demand_size 200 \
          --min_weight 0 \
          --max_weight 100 \
          --time 1000 \
          --beta $param \
          --trace_location traces/$trace &
      done
    elif [ $cache == 'LMIN' ]; then
      for param in .005 .01 .05 .1 .5 1 2 5; do
        echo python CacheNetwork.py res/$path/$path \
          --graph_type $graph \
          --cache_type $cache \
          --graph_size 100 \
          --min_capacity 1 \
          --max_capacity 5 \
          --query_nodes 5 \
          --catalog_size 50 \
          --demand_size 200 \
          --min_weight 0 \
          --max_weight 100 \
          --time 1000 \
          --gamma $param \
          --trace_location traces/$trace &
      done
    elif [ $cache == 'TBGRD' ]; then
      for freq in -1; do
        for param in .0001 .0005 .001 .005 .01 .05 .1; do
          echo python CacheNetwork.py res/$path/$path \
            --graph_type $graph \
            --cache_type $cache \
            --graph_size 100 \
            --min_capacity 1 \
            --max_capacity 5 \
            --query_nodes 5 \
            --catalog_size 50 \
            --demand_size 200 \
            --min_weight 0 \
            --max_weight 100 \
            --time 1000 \
            --action_selector_eta $param \
            --frequency $freq \
            --trace_location traces/$trace &
        done
      done
    else
      echo python CacheNetwork.py res/$path/$path \
        --graph_type $graph \
        --cache_type $cache \
        --graph_size 100 \
        --min_capacity 1 \
        --max_capacity 5 \
        --query_nodes 5 \
        --catalog_size 50 \
        --demand_size 200 \
        --min_weight 0 \
        --max_weight 100 \
        --time 1000 \
        --trace_location traces/$trace &
    fi
  done
done
sleep 1000