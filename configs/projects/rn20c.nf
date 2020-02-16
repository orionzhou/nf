env.NXF_HOME = '/home/springer/zhoux379/projects/nf/run/rn20c'

includeConfig '/home/springer/zhoux379/projects/nf/configs/nextflow.config'
includeConfig '/home/springer/zhoux379/projects/nf/configs/rnaseq.config'

params {
  profile = 'msi'
  genome = 'Zmays_B73'
  name = 'rn20c'
  input = "../../design/${params.name}.csv"
  singleEnd = false
  outdir = "../../raw/${params.name}"
}

process {
}

// Function to ensure that resource requirements don't go beyond
// a maximum limit
def check_max(obj, type) {
  if (type == 'memory') {
    try {
      if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
        return params.max_memory as nextflow.util.MemoryUnit
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'time') {
    try {
      if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
        return params.max_time as nextflow.util.Duration
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
      return obj
    }
  } else if (type == 'cpus') {
    try {
      return Math.min( obj, params.max_cpus as int )
    } catch (all) {
      println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
      return obj
    }
  }
}