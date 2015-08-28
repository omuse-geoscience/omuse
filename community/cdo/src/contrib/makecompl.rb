#!/usr/bin/env ruby
require 'optparse'

$opts = {:bin => '../src/cdo',:outfile => 'cdoCompletion'}
OptionParser.new do |o|
  o.banner = "Create auto completion files for different shells\n    " +
             "Usage:  makecompl.rb [-h] [-b CMD] [-o FILENAME]"
  o.separator("")
  o.separator("Options are ...")
  o.on("-b",
       "--binary CMD",
       "Choose a CDO binary different from #{$opts[:bin]}") {|cmd| $opts[:bin] = cmd}
  o.on("-o",
       "--outfile FILENAME",
       "Set the filename for all completions files (default: #{$opts[:outfile]})") {|fname| $opts[:outfile] = fname}
  o.on("-h","--help","Show this help") {puts o
    puts
    puts <<-'END'
RESTRICTIONS:
Supported shells are TCSH, BASH and ZSH. For tcsh completion is only performed
for regular options and operators with prepended '-'. Bash and zsh also
complete opertors without leading '-'.

AUTHOR:  Ralf Mueller, ralf.mueller@zmaw.de

LICENSE: CDO's License
END
    exit
  }
end.parse!
#=============================================================================== 
def getOperators
  # try to run the CDO binary first
  cmd       = $opts[:bin] + ' 2>&1'
  help      = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}
  if 5 >= help.size
    puts "Operators could not get listed by running the CDO binary (#{$opts[:bin]})"
    puts "Create operator list by scanning the documentation..."
    Dir.glob('../doc/tex/mod/*').map {|mod|
      File.open(mod).readlines.grep(/Operators/).map {|line| 
        line.chomp.split('=')[-1].split(' ')
      }
    }.flatten
  else
    help[(help.index("Operators:")+1)..help.index(help.find {|v| v =~ /CDO version/}) - 2].join(' ').split
  end
end
def getOptions
  cmd     = $opts[:bin] + ' 2>&1'
  options = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}.find_all {|item| /^-/.match(item)}.map {|o| o[0,2]}
  if options.empty?
    puts "Commandline options could not get listed by running the CDO binary (#{$opts[:bin]})"
    puts "Go on processing operators only ..."
    return []
  end
  options
end
#=============================================================================== 
operators = getOperators.sort
options   = getOptions.sort
# require 'pp'; pp operators; pp options
# Create the configuration files
complCmds = { 
  :tcsh => ['set cdoCmpl = (\\','); complete cdo \'c/-/$cdoCmpl/\' \'n/*/f/\''],
  :zsh  => ['compctl -k "('    ,')" -f cdo'],
  :bash => ['complete -W "'    ,'" -f cdo']
}
[:bash, :zsh, :tcsh].each {|shell|
  # tcsh: Remove the prepended '-' from the cdo cmdline options. This will be
  #       added through tcsh's completion command
  # otherwise: Add operators WITH leading '-'
  completions = (:tcsh == shell ) ? options.map {|o| o[1..-1]} + operators : options + operators.map {|o| "#{o} -#{o}"}
  File.open($opts[:outfile] + '.' + shell.to_s,'w') {|f|
    f << complCmds[shell][0] << "\n"
    completions.each {|item| f << item << " \\" << "\n" }
    f << complCmds[shell][1] << "\n"
  }
}
