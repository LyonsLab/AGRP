
import java.util.*;
import java.io.*;

public class TestGetContigOutputAndScaffoldInputNoSubgenomeInfo{
	public static void main( String[] args) throws Exception{
	
		String inputInfo = args[0];
        
		FileReader fr = new FileReader(inputInfo);
		BufferedReader br = new BufferedReader(fr);
        String aline = br.readLine();
        String[] info = aline.split("\t");
        int numberOfGenomes =Integer.parseInt(info[1]);
        System.out.println("numberOfGenomes\t"+numberOfGenomes);
     //   aline = br.readLine();
    
        
        int[] allGenomeIndex = new int[numberOfGenomes];
        int[] ploidyNumber = new int[numberOfGenomes];
      
        int[][] weightTable = new int[numberOfGenomes][ploidyNumber[0]];
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
        
        GenomeInString[] leaveGenomes = new GenomeInString[numberOfGenomes];
        String leaveGenomesFile = "outputFiles/genomesInString";
        for(int i =0; i<  allGenomeIndex.length; i++){
            leaveGenomesFile = leaveGenomesFile+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        leaveGenomesFile = leaveGenomesFile+".txt";
        
        FileReader fr1 = new FileReader(leaveGenomesFile);
        BufferedReader br1 = new BufferedReader(fr1);
        String aline1 = br1.readLine();
        info = aline1.split("\t");
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
        String rangeFile ="outputFiles/subgenomeRangesInGeneOrder_8400_9050_10997_19515.txt";
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
		String inputFile = "outputFiles/contigOutput.txt";
		//System.out.println(inputFile);
		//int seed = 2;
		//System.out.println("seed to reorder gene content\t"+ seed);
		
		fr1 = new FileReader(inputFile);
		br1 = new BufferedReader(fr1);
		aline = br1.readLine();
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
        String outPutFileName = "outputFiles/contig";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        outPutFileName = outPutFileName+".txt";
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
        
        aline = br.readLine();
        info = aline.split("\t");
        int minimumLength = Integer.parseInt(info[1]);
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
        aline = br.readLine();
        info = aline.split("\t");
        int ancChrNumber = Integer.parseInt(info[1]);
		
        int totalPloidy = 0;
        for(int i = 0; i< ploidyNumber.length; i++){
            totalPloidy = totalPloidy+ploidyNumber[i];
        }
         int[] weightTable2 =new int[totalPloidy];
         aline = br.readLine();
        for(int i = 0; i< totalPloidy;i++){
            aline = br.readLine();
            System.out.println(aline);
            info = aline.split("\t");
            weightTable2[i] = Integer.parseInt(info[1]);
            System.out.println("weight ="+weightTable2[i]);
        }
        for(int i = 0; i< weightTable2.length; i++){
            System.out.print(weightTable2[i]+"\t");
        }
        System.out.println();
       
        
        //**************
        //**************
        
       // for(int colorCode = 1; colorCode<ancChrNumber+1; colorCode++){
            GenomeInString[] genomeInContig = contigData[0].getGenomesInContig(contigData, -1, rangesForGenomes);
            //printout
            outPutFileName = "outputFiles/leaveGenomesInContig.txt";
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
            
            outPutFileName = "outputFiles/scaffoldInput.py";
            fstream = new FileWriter(outPutFileName,false);
            fbw = new BufferedWriter(fstream);
           // fbw.write("\"\"\"\n");
            String fileName = "mwmatching.py";
            fbw.write("import time\n");
            fbw.write("import sys\n");
            fbw.write("start = time.clock()\n");
            fbw.write("sys.setrecursionlimit(4000)\n");
            fr = new FileReader(fileName);
            br = new BufferedReader(fr);
            aline = br.readLine();
            while(aline!=null){
                fbw.write(aline+"\n");
                aline = br.readLine();
            }
            fbw.write("maxWeightMatching([ ");
            for(int i = 0; i< edgeMatrix.length-1; i++){
                for(int j = i+1; j< edgeMatrix.length; j++){
                    if(edgeMatrix[i][j]!=0){
                        fbw.write("("+i+","+j+","+edgeMatrix[i][j]+"),");
                    }
                }
            }
            fbw.write("])\n");
            fbw.write("end = time.clock()\n");
            fbw.write("print end - start\n");
            fbw.flush();
            fbw.close();

            
            
        

        
    
	}
}