
import java.util.*;
import java.io.*;

public class TestScaffoldInput{
	public static void main( String[] args) throws Exception{
        
        // read in subgenome inforamtion
        int numberOfGenomes = 4;
        SubgenomeRanges[] rangesForGenomes = new SubgenomeRanges[numberOfGenomes];
        String rangeFile ="outputFiles/subgenomeRangesInGeneOrder_8400_9050_10997_19515.txt";
        FileReader fr2 = new FileReader(rangeFile);
        BufferedReader br2 = new BufferedReader(fr2);
        String aline2 = "";
        String[] info = new String[0];
        for(int g = 0; g< rangesForGenomes.length; g++){
            aline2 = br2.readLine();
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

        // read in contig
        String[] root2 = new String[2765];
        String contigFile = "outputFiles/Checkcontig_8400_9050_10997_19515.txt";
        fr2 = new FileReader(contigFile);
        br2 = new BufferedReader(fr2);
        for(int g = 0; g< root2.length; g++){
            aline2 = br2.readLine();
            root2[g] = br2.readLine();
        }
        // read in leaveGenomes
        String leaveGenomesFile = "outputFiles/genomesInString_8400_9050_10997_19515.txt";
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
        
        GenomeInString[] leaveGenomes = new GenomeInString[4];
        
        for(int i = 0; i< leaveGenomes.length; i++){
            aline1 = br1.readLine();  // genome name
            String[] ag = new String[chrNumber[i]];
            for(int j = 0; j< ag.length; j++){
                aline1 = br1.readLine();
                ag[j] = br1.readLine();
            }
            leaveGenomes[i] = new GenomeInString(ag);
        }

        
		// for each contig, get all positions and gene groups
        int[] ploidy4 = new int[4]; ploidy4[0] = 3; ploidy4[1] = 3; ploidy4[2] = 3; ploidy4[3] = 1;
        int[] allGenomeIndex = new int[4]; allGenomeIndex[0] = 8400;
        allGenomeIndex[1] = 9050;
        allGenomeIndex[2] = 10997;
        allGenomeIndex[3] = 19515;
		Contig[] contigFL = new Contig[root2.length];
		for(int i = 0; i< contigFL.length; i++){
       // for(int i = 561; i< 562; i++){
        
			contigFL[i] = new Contig(i, root2[i],ploidy4, leaveGenomes,allGenomeIndex);
            for(int j = 0; j< numberOfGenomes; j++){
                if(rangesForGenomes[j].ranges.length!=0){
                    contigFL[i].getGeneGroupInAGenome(3, j,rangesForGenomes[j].ranges);
                }else{
                    contigFL[i].getGeneGroupNoSubgenomeInfoInAGenome(1,j,rangesForGenomes[j].ranges);
                }
            }
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
     
		
		int colorCode = 7;
		System.out.println("  for color code "+ colorCode);
        
        GenomeInString[] genomeInContig = contigFL[0].getGenomesInContig(contigFL, colorCode, rangesForGenomes);
				
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