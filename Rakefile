require 'bundler/setup'
$:.unshift File.expand_path('../lib', __FILE__)

# rake spec
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) { |t| t.verbose = false   }

# rake console
task :console do
  require 'pry'
  require 'gem_name'
  ARGV.clear
  Pry.start
end

desc "test simulations by running just 1 job"
task :test do
  system('qR.rb -j simtest -r -p skylake -t "02:00:00" -c 6 -y 0-0 simulate.R')
end

desc "run simulations"
task :sim do
  system('qR.rb -j sim -r -p skylake-himem -t "02:00:00" -c 6 -y 1-999 simulate.R')
end

desc "collate results"
task :collate do
  system('qR.rb -j collate -r -p skylake -c 6 -t "02:00:00" collate.R')
end

desc "show missing"
task :missing do
  wanted= (0..999).to_a.map { |i| "data/output-final/out#{i}.csv.gz" }
  finalfiles=Dir.glob('data/output-final/*')
  missing = (wanted - finalfiles).map { |f| File.basename(f,".csv.gz").sub("out","").to_i() }
  missing.each{ |i|
    puts "SLURM_ARRAY_TASK_ID=#{i} Rscript simulate.R"
  }
  puts "qR.rb -r -c 6 -p skylake-himem -j simcatch -y #{missing.join(",")} simulate.R"
end


desc "describe status"
task :status do
  gemmafiles=Dir.glob('data/output-new/*E.assoc.txt.gz')
  finalfiles=Dir.glob('data/output-final/*')
  puts "GEMMA runs complete: #{gemmafiles.length}"
  puts "sim runs complete: #{finalfiles.length}"
end

require 'fileutils'

desc "clean gemma output"
task :clean do
  FileUtils.rm Dir.glob('data/output-new/*')
end

desc "clean everything, to enable a fresh run"
task :deepclean => :clean do
  FileUtils.rm Dir.glob('data/output-final/*')
end
