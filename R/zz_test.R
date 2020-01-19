`%AND%` <- function(map1, map2) map1 & map2
`%OR%` <- function(map1, map2) map1 | map2
`%XOR%` <- function(map1, map2) (map1 & !map2) | (!map1 & map2)
`%MASK%` <- function(map1, map2) map1 * map2
AND <- function(map1, map2) map1 %AND% map2
OR <- function(map1, map2) map1 %OR% map2
XOR <- function(map1, map2) map1 %XOR% map2
MASK <- function(map1, map2) map1 %MASK% map2