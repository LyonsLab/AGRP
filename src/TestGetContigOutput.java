
import java.util.*;
import java.io.*;

public class TestGetContigOutput{
	public static void main( String[] args) throws Exception{
		System.out.println("\"\"\"");
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
        for(int i = 0; i< contigs.length; i++){
            fbw.write("contig "+ i+"\n"+contigs[i]+"\n");
		}
        
        fbw.flush();
        fbw.close();
    
	}
}