
import java.util.*;
import java.io.*;

public class TestScaffoldInputOLdData{
	public static void main( String[] args) throws Exception{
        
        System.out.println("\"\"\"");
		System.out.println( "genome index: 0: amborella, 1:peach,2: cacao 3:vitis");
		
		FileReader fr1;
		BufferedReader br1;
		String fn = "dataForTesting/4genomesAPCV.txt";
		System.out.println(fn);
		
		fr1 = new FileReader(fn);
		br1 = new BufferedReader(fr1);
		
		String[] g1 = new String[128];   // amborella
		String[] g2 = new String[8];   // peach
		String[] g3 = new String[10];   // cacao
		String[] g4 = new String[19];    // vitis
		
		
		// get Gene Content
		String geneContent = "";
		for(int i = 1; i<12834; i++){
			geneContent = geneContent+"  "+ new Integer(i).toString();
		}
		
		String aline = br1.readLine();aline = br1.readLine();aline = br1.readLine();aline = br1.readLine();
		aline = br1.readLine();for(int i = 0; i< g1.length; i++){aline = br1.readLine();g1[i] = br1.readLine();}
		aline = br1.readLine();for(int i = 0; i< g2.length; i++){aline = br1.readLine();g2[i] = br1.readLine();}
		aline = br1.readLine();for(int i = 0; i< g3.length; i++){aline = br1.readLine();g3[i] = br1.readLine();}
		aline = br1.readLine();for(int i = 0; i< g4.length; i++){aline = br1.readLine();g4[i] = br1.readLine();}
		
		
		GenomeInString[] allgenomes= new GenomeInString[4];
		
		allgenomes[0] = new GenomeInString(g1);
		allgenomes[1] = new GenomeInString(g2);
		allgenomes[2] = new GenomeInString(g3);
		allgenomes[3] = new GenomeInString(g4);
		
		
        // readIn ranges in gene order
		
		int[][] rangesInGenome1 = new int[30][5];
		int[][] rangesInGenome2 = new int[30][5];
		int[][] rangesInGenome3 = new int[30][5];
		aline = br1.readLine();
		for(int i = 0; i< rangesInGenome1.length; i++){
			aline = br1.readLine();
			String[] vs = aline.split("\t");
			for(int j = 0; j< vs.length; j++){
				rangesInGenome1[i][j] = new Integer(vs[j]).intValue();
			}
			rangesInGenome1[i][2] = rangesInGenome1[i][2]-1;
		}
		aline = br1.readLine();
		for(int i = 0; i< rangesInGenome2.length; i++){
			aline = br1.readLine();
			String[] vs = aline.split("\t");
			for(int j = 0; j< vs.length; j++){
				rangesInGenome2[i][j] = new Integer(vs[j]).intValue();
			}
			rangesInGenome2[i][2] = rangesInGenome2[i][2]-1;
		}
		aline = br1.readLine();
		for(int i = 0; i< rangesInGenome3.length; i++){
			aline = br1.readLine();
			String[] vs = aline.split("\t");
			for(int j = 0; j< vs.length; j++){
				rangesInGenome3[i][j] = new Integer(vs[j]).intValue();
			}
			rangesInGenome3[i][2] = rangesInGenome3[i][2]-1;
		}
		
		for(int i = 0; i< rangesInGenome1.length; i++){
			System.out.println(rangesInGenome1[i][0]+"\t"+rangesInGenome1[i][1]+"\t"+rangesInGenome1[i][2]+"\t"+rangesInGenome1[i][3]+"\t"+rangesInGenome1[i][4]);
		}
		for(int i = 0; i< rangesInGenome2.length; i++){
			System.out.println(rangesInGenome2[i][0]+"\t"+rangesInGenome2[i][1]+"\t"+rangesInGenome2[i][2]+"\t"+rangesInGenome2[i][3]+"\t"+rangesInGenome2[i][4]);
		}
		for(int i = 0; i< rangesInGenome3.length; i++){
			System.out.println(rangesInGenome3[i][0]+"\t"+rangesInGenome3[i][1]+"\t"+rangesInGenome3[i][2]+"\t"+rangesInGenome3[i][3]+"\t"+rangesInGenome3[i][4]);
		}
        
		
		// read in adjacencies
		String[] edges = new String[30000];
		int index = 0;
		String inputFile = "dataForTesting/APCVFirstLevelContigAW160Adj.txt";
		System.out.println(inputFile);
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
		int weightToInclude = 1;
		System.out.println("edge with Weight >="+weightToInclude);
		GenomeAdj allEdges2 = new GenomeAdj(geneContent);
        //	allEdges2.reorderGeneContent(seed);
		//System.out.println(allEdges2.rootGeneContent);
		
		allEdges2.initialValue(edges);
		String[] root2 = allEdges2.getGenomeHighEdgeWeight(weightToInclude);
		System.out.println("root Genome (in genes )");
		for(int i = 0; i< root2.length; i++){
            //	System.out.println("chr "+ i+"\n"+root2[i]);
		}
        
		// for each contig, get all positions and gene groups
        int[] ploidy4 = new int[4]; ploidy4[0] = 1; ploidy4[1] = 3; ploidy4[2] = 3; ploidy4[3] = 3;
        int[] allGenomeIndex = new int[4]; allGenomeIndex[0] = 19515;
        allGenomeIndex[1] = 8400;
        allGenomeIndex[2] = 10997;
        allGenomeIndex[3] = 9050;
		Contig[] contigFL = new Contig[root2.length];
		for(int i = 0; i< contigFL.length; i++){
       // for(int i = 561; i< 562; i++){
            System.out.println("for contig" + root2[i]);
			contigFL[i] = new Contig(i, root2[i],ploidy4, allgenomes,allGenomeIndex);
           /* for(int j = 0; j< 4; j++){
                if(rangesForGenomes[j].ranges.length!=0){
                    contigFL[i].getGeneGroupInAGenome(3, j,rangesForGenomes[j].ranges);
                }else{
                    contigFL[i].getGeneGroupNoSubgenomeInfoInAGenome(1,j,rangesForGenomes[j].ranges);
                }
            }*/
            contigFL[i].getGeneGroupInAGenome(3,1,rangesInGenome1);
            contigFL[i].getGeneGroupInAGenome(3,2,rangesInGenome2);
            contigFL[i].getGeneGroupInAGenome(3,3,rangesInGenome3);
            
          //  contigFL[i].printAContig();
		}

        // count how many contigs with more than 1 colorCode, and how many contigs no geneGroup
		int totalContigNoGG = 0;
		int[] contigNoGGLength = new int[2000];
		System.out.println("======contigs without geneGroups");
		for(int i = 0; i< contigFL.length; i++){
			if(contigFL[i].ggIndex == 0){
				totalContigNoGG++;
				contigNoGGLength[contigFL[i].totalGeneNumber]++;
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
		for(int i = 0; i< contigFL.length; i++){
			boolean twoOrMore = false;
			for(int j = 0; j< contigFL[i].ggIndex-1; j++){
				for(int k= j+1; k<contigFL[i].ggIndex; k++){
					if(contigFL[i].geneGroups[j].colorCode != contigFL[i].geneGroups[k].colorCode){
						twoOrMore= true;
					}
				}
			}
			if(twoOrMore == true){
				totalContigMoreColor++;
				contigMoreColorLength[contigFL[i].totalGeneNumber]++;
				//contigFL[i].printAContig();
			}
		}
		
		System.out.println(" contigs with GeneGroups in more than one color  "+ totalContigMoreColor);
		System.out.println(" the length of contig with GeneGroups in more than one color ");
		for(int i = 0; i< contigMoreColorLength.length; i++){
			if(contigMoreColorLength[i]!=0){
				System.out.println(i+"\t"+ contigMoreColorLength[i]);}
		}
     
		
	//	int colorCode = 7;
	//	System.out.println("  for color code "+ colorCode);
        
    //    GenomeInString[] genomeInContig = contigFL[0].getGenomesInContig(contigFL, colorCode, rangesForGenomes);
				
		/*int weightAmborella = 100;
		int weightPCV = 150;
		int weightCVorPV = 125;
		int weightPC = 110;
		int weightPorCorV = 100;
		System.out.println("weightAmborella "+ weightAmborella);
		
		EdgeForMWM cc1Edge = new EdgeForMWM(genomeInContig, genomeInContig[0].chrs[0]);
		byte[][] edgeInAmborella = cc1Edge.getMWMInput(genomeInContig[1].chrs,1);
		
		byte[][] edgeInPeach1 = cc1Edge.getMWMInput(genomeInContig[2].chrs,1);
		byte[][] edgeInCacao1 = cc1Edge.getMWMInput(genomeInContig[3].chrs,1);
		byte[][] edgeInVitis1 = cc1Edge.getMWMInput(genomeInContig[4].chrs,1);
		
		byte[][] edgeInPeach2 = cc1Edge.getMWMInput(genomeInContig[5].chrs,1);
		byte[][] edgeInCacao2 = cc1Edge.getMWMInput(genomeInContig[6].chrs,1);
		byte[][] edgeInVitis2 = cc1Edge.getMWMInput(genomeInContig[7].chrs,1);
		
		byte[][] edgeInPeach3 = cc1Edge.getMWMInput(genomeInContig[8].chrs,1);
		byte[][] edgeInCacao3 = cc1Edge.getMWMInput(genomeInContig[9].chrs,1);
		byte[][] edgeInVitis3 = cc1Edge.getMWMInput(genomeInContig[10].chrs,1);
		
		
		int[][] edgeMatrix = new int[edgeInAmborella.length][edgeInAmborella.length];
		for(int i = 0; i< edgeInAmborella.length; i++){
			for(int j = 0; j< edgeInAmborella.length; j++){
				edgeMatrix[i][j] = edgeMatrix[i][j]+edgeInAmborella[i][j]*weightAmborella;
				//light
				if(edgeInPeach1[i][j]==1 && edgeInCacao1[i][j]==1 && edgeInVitis1[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPCV;}
	
				if(edgeInPeach1[i][j]==0 && edgeInCacao1[i][j]==1 && edgeInVitis1[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				if(edgeInPeach1[i][j]==1 && edgeInCacao1[i][j]==0 && edgeInVitis1[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				
				if(edgeInPeach1[i][j]==1 && edgeInCacao1[i][j]==1 && edgeInVitis1[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPC;}
				
				if(edgeInPeach1[i][j]==0 && edgeInCacao1[i][j]==0 && edgeInVitis1[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach1[i][j]==0 && edgeInCacao1[i][j]==1 && edgeInVitis1[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach1[i][j]==1 && edgeInCacao1[i][j]==0 && edgeInVitis1[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				//=============medium
				if(edgeInPeach2[i][j]==1 && edgeInCacao2[i][j]==1 && edgeInVitis2[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPCV;}
				
				if(edgeInPeach2[i][j]==0 && edgeInCacao2[i][j]==1 && edgeInVitis2[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				if(edgeInPeach2[i][j]==1 && edgeInCacao2[i][j]==0 && edgeInVitis2[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				
				if(edgeInPeach2[i][j]==1 && edgeInCacao2[i][j]==1 && edgeInVitis2[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPC;}
				
				if(edgeInPeach2[i][j]==0 && edgeInCacao2[i][j]==0 && edgeInVitis2[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach2[i][j]==0 && edgeInCacao2[i][j]==1 && edgeInVitis2[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach2[i][j]==1 && edgeInCacao2[i][j]==0 && edgeInVitis2[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				
				//=============dark
				if(edgeInPeach3[i][j]==1 && edgeInCacao3[i][j]==1 && edgeInVitis3[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPCV;}
				
				if(edgeInPeach3[i][j]==0 && edgeInCacao3[i][j]==1 && edgeInVitis3[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				if(edgeInPeach3[i][j]==1 && edgeInCacao3[i][j]==0 && edgeInVitis3[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightCVorPV;}
				
				if(edgeInPeach3[i][j]==1 && edgeInCacao3[i][j]==1 && edgeInVitis3[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPC;}
				
				if(edgeInPeach3[i][j]==0 && edgeInCacao3[i][j]==0 && edgeInVitis3[i][j]==1){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach3[i][j]==0 && edgeInCacao3[i][j]==1 && edgeInVitis3[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				if(edgeInPeach3[i][j]==1 && edgeInCacao3[i][j]==0 && edgeInVitis3[i][j]==0){edgeMatrix[i][j] = edgeMatrix[i][j]+weightPorCorV;}
				
				
				
				
				
				
			}
		}
	
		*/
	/*	FileReader fr;
		BufferedReader br;
		System.out.println("\"\"\"");
		String fileName = "mwmatching.py";
		System.out.println("import time");
		System.out.println("import sys");
		
		System.out.println("start = time.clock()");
		System.out.println("sys.setrecursionlimit(4000)");
		
		fr = new FileReader(fileName);
		br = new BufferedReader(fr);
		aline = br.readLine();
		while(aline!=null){
			System.out.println(aline);
			aline = br.readLine();
		}
		System.out.print("maxWeightMatching([ ");
		for(int i = 0; i< edgeMatrix.length-1; i++){
			for(int j = i+1; j< edgeMatrix.length; j++){
				if(edgeMatrix[i][j]!=0){
					System.out.print("("+i+","+j+","+edgeMatrix[i][j]+"),");
				}
			}
		}
		System.out.println("])");
		
		System.out.println("end = time.clock()");
		
		System.out.println("print end - start");
		//time = Calendar.getInstance();
		//System.out.println("end " + time.getTime());
	*/
        
	}
}