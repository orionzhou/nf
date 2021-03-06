env.NXF_HOME = '/home/springer/zhoux379/projects/nf/run/rn20d2'

includeConfig '/home/springer/zhoux379/projects/nf/configs/nextflow.config'
// includeConfig '/home/springer/zhoux379/projects/nf/configs/smrnaseq.config'

params {
  profile = 'msi'
  genome = 'Zmays_B73'
  name = 'rn20d2'
  input = "../../design/${params.name}.csv"
  outdir = "../../raw/${params.name}"

  protocol = 'illumina'
  mature = '/home/springer/zhoux379/projects/nf/misc/mirbase/mature.fa.gz'
  hairpin = '/home/springer/zhoux379/projects/nf/misc/mirbase/hairpin.fa.gz'
  mirtrace_species = 'zma'
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