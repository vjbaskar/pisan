
# BASH based software suite
This is a wrap around for many of commonly used pipelines and supports submission through IBM LSF bsub command.
If you want to add other batch submission systems, then you have edit the <installdir>/bin/farmsub.sh file to support it.

# Installation Instructions
This is primarily written in linux bash. Some of these may be incompatible with Mac bash and needless to say Windows is out of question.
It requires you to install all the packages (some are given below) that it uses in its various sub commands.


The main config file in <installdir> ~/$HOME/.pisan/conf
If you provide a conf file with the commands it will overwrite options of pisan.conf

In the conf file you have to change the variable to point to the path where pisan is installed 
`mconf_installdir="$HOME/Software/pisan"`

```
echo "PATH=$USER/Software/pisan/bin:$PATH" >> $HOME/.bashrc
```

# Main command

`pisan` - Generates all the suite commands available
each command is a bash file in the bin folder.


# Depends on
All relevant tools should be added to your $PATH. The main ones are

* Rscript
* BASH
* BASH -> getopts
* samtools
* biobambam2
* anaconda python = create a python2.7 environment with name p27
* anaconda python = create a python3 environment with name p3



# How it works

pisan sources the main and local config files. Local config files overwrite main config

It then reads in arguments from src/general_args.sh and functions from src/basic_functions.sh

