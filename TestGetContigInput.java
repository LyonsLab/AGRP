
import java.util.*;
import java.io.*;

public class TestGetContigInput{
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
       // aline = br.readLine();
      //  System.out.println("genomeFies\t"+aline);
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
        
      /*  for(int i = 0; i< leaveGenomes.length; i++){
            System.out.println("a genome \t"+allGenomeIndex[i]);
            leaveGenomes[i].print();
        }*/
		// get Gene Content
		String geneContent = "";
		for(int i = 1; i<numberOfGenes+1; i++){
			geneContent = geneContent+"  "+ new Integer(i).toString();
		}
	
		EdgeForMWM allEdges = new EdgeForMWM(leaveGenomes, geneContent);
        
        int[][] edgeMatrix = new int[allEdges.nodeString.length][allEdges.nodeString.length];
        
        for(int i = 0; i< leaveGenomes.length; i++){
            byte[][] edgeMatrix1 = allEdges.getASetOfEdge(leaveGenomes[i].chrs);
            for(int j = 0; j< edgeMatrix.length; j++){
                for(int k = 0; k< edgeMatrix.length; k++){
                    int weighthere = weightTable[i][edgeMatrix1[j][k]];
                    edgeMatrix[j][k] = edgeMatrix[j][k] +weighthere;
                }
            }
        }
        
       /* for(int i = 0; i< edgeMatrix.length-1; i++){
            for(int j = i+1; j< edgeMatrix.length; j++){
                if(edgeMatrix[i][j]>0){
                    System.out.println(i+"\t"+j+"\t"+allEdges.nodeString[i]+"\t"+allEdges.nodeString[j]+"\t"+edgeMatrix[i][j]);
                }
            }
        }
     */
        
        /*int[] edgeWeight = new int[2000];
		for(int i = 0; i< edgeMatrix.length-1; i++){
			for(int j = i+1; j< edgeMatrix.length; j++){
				edgeWeight[edgeMatrix[i][j]]++;
			}
		}
	
		System.out.println("Weight\tnumberOFAdj");
		for(int i = 0; i< edgeWeight.length; i++){
			if(edgeWeight[i]!=0){
				System.out.println(i+ "\t" + edgeWeight[i]);
			}
		}
		
		int totalNumberOfEdge = 0;
		for(int i = 0; i< edgeMatrix.length-1; i++){
			for(int j = i+1; j< edgeMatrix.length; j++){
				if(edgeMatrix[i][j]!=0){
					totalNumberOfEdge++;
				}
			}
		}*/
        
        
        String outPutFileName = "outputFiles/contigInput";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
      /*  outPutFileName = outPutFileName+".py";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
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
                if(edgeMatrix[i][j]>=weightToInclude){
                    fbw.write("("+i+","+j+","+edgeMatrix[i][j]+"),");
                }
            }
        }
        fbw.write("])\n");
        fbw.write("end = time.clock()\n");
        fbw.write("print end - start\n");
       */
        outPutFileName = outPutFileName+".txt";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);
        for(int i = 0; i< edgeMatrix.length-1; i++){
            for(int j = i+1; j< edgeMatrix.length; j++){
                if(edgeMatrix[i][j]>=weightToInclude){
                    fbw.write(i+"\t"+j+"\t"+edgeMatrix[i][j]+"\n");
                }
            }
        }
        fbw.flush();
        fbw.close();
    
    
	}
}