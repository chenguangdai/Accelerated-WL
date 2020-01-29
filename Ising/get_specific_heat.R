### calculate the specific heat
get_specific_heat <- function(logdensity, energy, temperature){
  logw <- logdensity - energy/temperature
  maxlw <- max(logw)
  logw <- logw - maxlw
  w <- exp(logw)/sum(exp(logw))
  average_energy <- sum(energy * w)
  average_squared_energy <- sum(energy^2 * w)
  specific_heat <- (average_squared_energy - average_energy^2)/temperature^2
  return(specific_heat)
}
