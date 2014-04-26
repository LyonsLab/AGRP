
import java.util.*;
import java.io.*;
import java.util.regex.Pattern;

public class TestGetContigOutputAndScaffoldInput{
	public static void main( String[] args) throws Exception{
        if (args.length < 8) {
            System.out.println(Usage());
            System.exit(1);
        }
        
        String[] genomeIndex = null;
        String[] weights = null;
        String[] ploidyArray = null;
        int weightToInclude = 1;
        String minimumLengthString = null;
        String inputGenomes = null;
        String output_directory = ".";
 
        // Input files
        String contigOutputFile = null;
        String rangeFile = null;
        String genomeInStringFile = null;
        
        for(int i = 0; i < args.length; i += 2) {
            String argument = args[i].substring(1);
            String option = args[i + 1];
            
            switch(argument) {
                case "mml":
                    minimumLengthString = option;
                case "p":
                    String ploidyString = option;
                    ploidyArray = ploidyString.split(Pattern.quote(","));
                    break;
                case "w":
                    String weightSchemeString = option;
                    weights = option.split(Pattern.quote(","));
                    break;
                case "g":
                    String genomeIndexString = option;
                    genomeIndex = option.split(Pattern.quote(","));
                    break;
                case "wa":
                    weightToInclude = Integer.parseInt(option);
                    break;
                case "co":
                    File contig_file = new File(option);
                    if (!contig_file.exists()) {
                        System.out.println("The contig output directory does not exist.");
                        System.exit(1);
                    }
                    contigOutputFile = option;
                    break;
                case "s":
                    File sub_genome_file = new File(option);
                    if (!sub_genome_file.exists()) {
                        System.out.println("The sub genome file does not exist.");
                        System.exit(1);
                    }
                    rangeFile = option;
                    break;
                case "gf":
                    File genome_file = new File(option);
                    if (!genome_file.exists()) {
                        System.out.println("The genomeInString file does not exist.");
                        System.exit(1);
                    }
                    genomeInStringFile = option;
                    break;
                case "o":
                    output_directory = option;
                    break;
                default:
                    System.out.println("Unknown argument");
                    System.out.println(Usage());
                    System.exit(1);
            }
        }
        
        if (genomeIndex == null || weights == null || minimumLengthString == null ||
            contigOutputFile == null || rangeFile == null || ploidyArray == null ||
            genomeInStringFile == null) {
            System.out.println("Please specify all the command line arguments");
            System.out.println(Usage());
            System.exit(1);
        }
        
        File fd1 = new File(output_directory);
        if (!fd1.exists()) {
            fd1.mkdir();
        }
        
        File fd2 = new File(output_directory + "/scaffolds");
        if (!fd2.exists()) {
            fd2.mkdir();
        }
        
        File fd3 = new File(output_directory + "/binfiles");
        if (!fd3.exists()) {
            fd3.mkdir();
        }
        
	/*	FileReader fr = new FileReader(inputInfo);
		BufferedReader br = new BufferedReader(fr);
        String aline = br.readLine();
        String[] info = aline.split("\t");
        int numberOfGenomes =Integer.parseInt(info[1]);
        System.out.println("numberOfGenomes\t"+numberOfGenomes);
     //   aline = br.readLine();
    */
        
        
        int numberOfGenomes = genomeIndex.length;
        int[] allGenomeIndex = new int[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            allGenomeIndex[i] = Integer.parseInt(genomeIndex[i]);
        }
        
        int[] ploidyNumber = new int[numberOfGenomes];  ///**************
        for(int i = 0; i< ploidyArray.length; i++){
            ploidyNumber[i] =Integer.parseInt(ploidyArray[i]);
        }
     /*   int[][] weightTable = new int[numberOfGenomes][ploidyNumber[0]];
        aline = br.readLine();//genomeIndex	ploidyNumber	chrNumber	weights
        for(int i = 0; i< numberOfGenomes; i++){
            aline = br.readLine();
            String[] infos = aline.split("\t");
            allGenomeIndex[i] = Integer.parseInt(infos[0]);
            ploidyNumber[i] = Integer.parseInt(infos[1]);
          //  chrNumber[i] = Integer.parseInt(infos[2]);
            weightTable[i] = new int[ploidyNumber[i]+1];
            for(int j = 0; j< weightTable[i].length; j++){
                weightTable[i][j] = Integer.parseInt(infos[2+j]);
            }
        }
      
        for(int i = 0; i< allGenomeIndex.length; i++){
            System.out.print(allGenomeIndex[i]+"\t"+ploidyNumber[i]+"\t");
            for(int j = 0; j< weightTable[i].length; j++ ){
                System.out.print(weightTable[i][j]+"\t");
            }
            System.out.println();
        }
        aline = br.readLine();
        info = aline.split("\t");
        int weightToInclude = Integer.parseInt(info[1]);
		System.out.println("edge with Weight >="+weightToInclude);
        */
        GenomeInString[] leaveGenomes = new GenomeInString[numberOfGenomes];
        String leaveGenomesFile = genomeInStringFile;
       
        FileReader fr1 = new FileReader(leaveGenomesFile);
        BufferedReader br1 = new BufferedReader(fr1);
        String aline1 = br1.readLine();
        String[] info = aline1.split("\t");
        int numberOfGenes =Integer.parseInt(info[1]);
        System.out.println("numberOfGenes\t"+numberOfGenes);
        int[] chrNumber = new int[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            aline1 = br1.readLine();
            System.out.println(aline1);
            info = aline1.split("\t");
            chrNumber[i]  =Integer.parseInt(info[3]);
        }
        
        for(int i = 0; i< leaveGenomes.length; i++){
            aline1 = br1.readLine();  // genome name
            String[] ag = new String[chrNumber[i]];
            for(int j = 0; j< ag.length; j++){
                aline1 = br1.readLine();
                ag[j] = br1.readLine();
            }
            leaveGenomes[i] = new GenomeInString(ag);
        }
        //
        SubgenomeRanges[] rangesForGenomes = new SubgenomeRanges[numberOfGenomes];
        
        FileReader fr2 = new FileReader(rangeFile);
        BufferedReader br2 = new BufferedReader(fr2);
        for(int g = 0; g< rangesForGenomes.length; g++){
           String aline2 = br2.readLine();
            info = aline2.split("\t");
            int gi = Integer.parseInt(info[0]);
            int b = Integer.parseInt(info[1]);
            int p = Integer.parseInt(info[2]);
            System.out.println(gi +"\t"+b+"\t"+p);
            int[][] r = new int[b][5];
            aline2 = br2.readLine();
            for(int i = 0; i< r.length; i++){
                aline2 = br2.readLine();
                info = aline2.split("\t");
                for(int j = 0; j< info.length; j++){
                    r[i][j] = Integer.parseInt(info[j]);
                }
            }
            rangesForGenomes[g] =new SubgenomeRanges(gi, p,b, r);
        }
        

 
		// get Gene Content
		String geneContent = "";
		for(int i = 1; i<numberOfGenes+1; i++){
			geneContent = geneContent+"  "+ new Integer(i).toString();
		}
        
        String[] edges = new String[numberOfGenes];
		int index = 0;
		
		//System.out.println(inputFile);
		//int seed = 2;
		//System.out.println("seed to reorder gene content\t"+ seed);
		
		fr1 = new FileReader(contigOutputFile);
		br1 = new BufferedReader(fr1);
		String aline = br1.readLine();
		while(aline!=null && aline.substring(0,5).equals("total")==false){
			edges[index] = aline;
			index++;
			aline = br1.readLine();
		}
		System.out.println("total edges "+ index);
		System.out.println("---------"+aline);
		
		// get genome (first level contigs)
		
		GenomeAdj allEdges2 = new GenomeAdj(geneContent);
        //	allEdges2.reorderGeneContent(seed);
		//System.out.println(allEdges2.rootGeneContent);
		
		allEdges2.initialValue(edges);
		String[] contigs = allEdges2.getGenomeHighEdgeWeight(weightToInclude);
        String outPutFileName = output_directory + "/contig2genes.txt";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        fbw.write("totalContig\t"+contigs.length+"\n");
        for(int i = 0; i< contigs.length; i++){
            fbw.write("contig "+ i+"\n"+contigs[i]+"\n");
		}
        
        fbw.flush();
        fbw.close();
        //*****************
        //*****************
        
      /*  aline = br.readLine();
        info = aline.split("\t");*/
        int minimumLength = Integer.parseInt(minimumLengthString);
		System.out.println("minimumLength >="+minimumLength);
        
        //*****************
        //*****************
        
        Contig[] contigData = new Contig[contigs.length];
        for(int i = 0; i< contigData.length; i++){
			contigData[i] = new Contig(i, contigs[i], ploidyNumber,leaveGenomes, allGenomeIndex);
			for(int j = 0; j< numberOfGenomes; j++){
                if(rangesForGenomes[j].ranges.length!=0){
                    contigData[i].getGeneGroupInAGenome(minimumLength,j,rangesForGenomes[j].ranges);
                }else{
                    contigData[i].getGeneGroupNoSubgenomeInfoInAGenome(minimumLength, j,rangesForGenomes[j].ranges);
                }
            }
           // contigData[i].printAContig();
		}
        
        // count how many contigs with more than 1 colorCode, and how many contigs no geneGroup
		int totalContigNoGG = 0;
		int[] contigNoGGLength = new int[2000];
		System.out.println("======contigs without geneGroups");
		for(int i = 0; i< contigData.length; i++){
			if(contigData[i].ggIndex == 0){
				totalContigNoGG++;
				contigNoGGLength[contigData[i].totalGeneNumber]++;
				//contigFL[i].printAContig();
			}
		}
		
		System.out.println(" contigs without GeneGroup   "+ totalContigNoGG);
		System.out.println(" the length of contigs without GeneGroup ");
		for(int i = 0; i< contigNoGGLength.length; i++){
			if(contigNoGGLength[i]!=0){
				System.out.println(i+"\t"+ contigNoGGLength[i]);}
		}
		int totalContigMoreColor = 0;
		int[] contigMoreColorLength = new int[2000];
		System.out.println("======contigs with geneGroups in more than one color");
		for(int i = 0; i< contigData.length; i++){
			boolean twoOrMore = false;
			for(int j = 0; j< contigData[i].ggIndex-1; j++){
				for(int k= j+1; k<contigData[i].ggIndex; k++){
					if(contigData[i].geneGroups[j].colorCode != contigData[i].geneGroups[k].colorCode){
						twoOrMore= true;
					}
				}
			}
			if(twoOrMore == true){
				totalContigMoreColor++;
				contigMoreColorLength[contigData[i].totalGeneNumber]++;
				contigData[i].printAContig();
			}
		}
		
		System.out.println(" contigs with GeneGroups in more than one color  "+ totalContigMoreColor);
		System.out.println(" the length of contig with GeneGroups in more than one color ");
		for(int i = 0; i< contigMoreColorLength.length; i++){
			if(contigMoreColorLength[i]!=0){
				System.out.println(i+"\t"+ contigMoreColorLength[i]);}
		}
     
        
        //****************
        //****************
     //   aline = br.readLine();
     //   info = aline.split("\t");
       // int ancChrNumber = Integer.parseInt(info[1]);
       
        int ancChrNumber = getAncChrNumber(rangesForGenomes);//*************
		
        int totalPloidy = 0;
        for(int i = 0; i< ploidyNumber.length; i++){
            totalPloidy = totalPloidy+ploidyNumber[i];
        }
        
        int[] weightTable2 =new int[totalPloidy];
        int wtIndex = 0;
        for(int i = 0; i< genomeIndex.length;i++){
            for(int j = 0; j< ploidyNumber[i]; j++){
                weightTable2[wtIndex] = Integer.parseInt(weights[i]);
                wtIndex++;
            }
        }
        for(int i = 0; i< weightTable2.length; i++){
            System.out.print(weightTable2[i]+"\t");
        }
        System.out.println();
       
        
        //**************
        //**************
        
        for(int colorCode = 1; colorCode<ancChrNumber+1; colorCode++){
            GenomeInString[] genomeInContig = contigData[0].getGenomesInContig(contigData, colorCode, rangesForGenomes);
            //printout
            outPutFileName = output_directory + "/binfiles/leaveGenomesInContigForAAncChr"+new Integer(colorCode)+".txt";
            fstream = new FileWriter(outPutFileName,false);
            fbw = new BufferedWriter(fstream);
            for(int ac= 0; ac<genomeInContig.length; ac++){
                fbw.write("genome\t"+ac+"\t"+genomeInContig[ac].chrs.length+"\n");
                for(int chrl = 0; chrl< genomeInContig[ac].chrs.length; chrl++){
                    fbw.write("chr\t"+chrl+"\n"+genomeInContig[ac].chrs[chrl]+"\n");
                }
            }
            fbw.flush();
            fbw.close();

            
         /*   String geneContentInContig = "";
            for(int i = 1; i<contigData.length; i++){
                geneContentInContig = geneContentInContig+"  contig"+ new Integer(i).toString();
            }*/
            
            EdgeForMWM allEdges = new EdgeForMWM(genomeInContig,genomeInContig[0].chrs[0]);
            
            int[][] edgeMatrix = new int[allEdges.nodeString.length][allEdges.nodeString.length];
            
            for(int i = 1; i< genomeInContig.length; i++){
                byte[][] edgeMatrix1 = allEdges.getASetOfEdge(genomeInContig[i].chrs);
                for(int j = 0; j< edgeMatrix.length; j++){
                    for(int k = 0; k< edgeMatrix.length; k++){
                        int weighthere = edgeMatrix1[j][k]*weightTable2[i-1];
                        edgeMatrix[j][k] = edgeMatrix[j][k] +weighthere;
                    }
                }
            }
            
            outPutFileName = output_directory + "/scaffolds/scaffoldInput"+new Integer(colorCode)+".txt";
            fstream = new FileWriter(outPutFileName,false);
            fbw = new BufferedWriter(fstream);
           // fbw.write("\"\"\"\n");
            for(int i = 0; i< edgeMatrix.length-1; i++){
                for(int j = i+1; j< edgeMatrix.length; j++){
                    if(edgeMatrix[i][j]!=0){
                        fbw.write(i+"\t"+j+"\t"+edgeMatrix[i][j]+"\n");
                    }
                }
            }
            fbw.flush();
            fbw.close();
        }
        
    }
    
   
       // int ancChrNumber = atest.getAncChrNumber(rangesForGenomes);//*************
        public static int getAncChrNumber(SubgenomeRanges[] allSubGenomeInfo){
            int result =0;
            for(int i = 0;i < allSubGenomeInfo.length; i++){
                for(int j = 0; j< allSubGenomeInfo[i].ranges.length; j++){
                    if(allSubGenomeInfo[i].ranges[j][0] > result){
                        result =allSubGenomeInfo[i].ranges[j][0];
                    }
                }
                
            }
            return result;
        }

        public static String Usage() {
            return "TODO";
        }
    
	
}