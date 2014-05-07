AGRP
----

This program is developed by Chunfeng Zheng, and refactored during Hackathon in
Tucson in April 2014. A more detailed description of discrete steps of this
pipeline is documented [here](http://genomevolution.org/wiki/index.php/Ancestral_Reconstruction_Pipeline).

##Compiling##
The source files can be built using ant. The command below will compile the files into the __build__ directory.

```bash
ant
```

##Run the pipeline##
Please check the pipeline in `run.sh`.

```bash
./run.sh
```

##Clean up##
To remove the build directory run:

```bash
ant clean
```
