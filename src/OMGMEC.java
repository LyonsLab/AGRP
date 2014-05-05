import java.util.*;
import java.io.*;

public class OMGMEC{
    Homolog[] originalHomologs;
    Homolog[] orthologSets;
    int[] allGenomeIndex;
    int genomeNumber;
    
    public OMGMEC(Homolog[] before, int[] allGenomeIndexHere){
        originalHomologs = new Homolog[before.length];
        originalHomologs = before;
        genomeNumber = allGenomeIndexHere.length;
        allGenomeIndex = new int[allGenomeIndexHere.length];
        allGenomeIndex = allGenomeIndexHere;
        
    }
    
    public void MECmethod(int biggraph, int lowSimilarity){
        int originalGood = 0;
        int maximum =0; for(int i = 0; i< originalHomologs.length; i++){maximum = maximum+originalHomologs[i].genePairs.length;}
        Homolog[] tmpOrtholog = new Homolog[maximum];
        int indexAll = 0;
        for(int i =0; i< originalHomologs.length; i++){
            Homolog ageneFamily = originalHomologs[i];
            if(ageneFamily.genePairs.length < biggraph){
                Homolog[] goodH = checkGoodOrNot(ageneFamily, allGenomeIndex); //*****************
				Homolog[] orthologSet = null;
				if(goodH.length==1){
					originalGood++;
					orthologSet = new Homolog[1];
					orthologSet=goodH;
                 
				}else{
					orthologSet= minimumDeletion(ageneFamily, allGenomeIndex);
				}
				for(int j = 0; j< orthologSet.length; j++){
					tmpOrtholog[indexAll] = orthologSet[j];
					indexAll++;
				}
                
            }else{
				GenePair[] edges = removeEdge(ageneFamily,lowSimilarity);
                
				Homolog[] allh = ageneFamily.getAllHomologs1(edges,genomeNumber);
				for(int j = 0; j< allh.length; j++){
					if(allh[j].genePairs.length<biggraph){
                        Homolog agenefamily2 = allh[j];
                        Homolog[] goodH2 = checkGoodOrNot(agenefamily2, allGenomeIndex);
                        Homolog[] orthologSet2 = null;
                        if(goodH2.length==1){
                            orthologSet2 = new Homolog[1];
                            orthologSet2 = goodH2;
                        }else{
                            orthologSet2= minimumDeletion(agenefamily2, allGenomeIndex);
                        }
                        for(int o = 0; o< orthologSet2.length; o++){
                            tmpOrtholog[indexAll] = orthologSet2[o];
                            indexAll++;
                        }
					}else{
						System.out.println("drop a graph with "+ allh[j].genePairs.length);
					}
				}
            }
			
        }
        System.out.println("original good"+originalGood);
        int[] homologSize1 = new int[10000];
        int[] homologSize2 = new int[10000];
        orthologSets = new Homolog[indexAll];
        for(int i = 0; i< orthologSets.length; i++){
            orthologSets[i] = tmpOrtholog[i];
        }
        
        for(int i = 0; i< orthologSets.length; i++){
            homologSize1[orthologSets[i].genePairs.length]++;
            //if(orthologSets[i].genePairs.length==38){
             //     orthologSets[i].print();
           // }
        }
        
        for(int i = 0; i< orthologSets.length; i++){
            homologSize2[orthologSets[i].genes.length]++;
        }
        System.out.println("size in gene pairs");
        for(int i = 0; i< homologSize1.length; i++){
            if(homologSize1[i]!=0){
                System.out.println(i+"\t"+homologSize1[i]);
            }
        }
        
        
        System.out.println("size in genes");
        for(int i = 0; i< homologSize2.length; i++){
            if(homologSize2[i]!=0){
                System.out.println(i+"\t"+homologSize2[i]);
            }
        }
        
        for(int i = 0; i< orthologSets.length-1; i++){
            for(int j = i+1; j< orthologSets.length; j++){
                for(int k = 0; k< orthologSets[i].genes.length; k++){
                    for(int h = 0; h< orthologSets[j].genes.length; h++){
                        if(orthologSets[i].genes[k].name.equals(orthologSets[j].genes[h].name)){
                            System.out.println("more combine "+i+"\t"+j);
                        }
                    }
                }
            }
        }
        
        // order them in size
        for(int i = 0; i< orthologSets.length-1; i++){
            for(int j = i+1; j<orthologSets.length; j++ ){
                if(orthologSets[i].genePairs.length< orthologSets[j].genePairs.length){
                    Homolog atmp = orthologSets[i];
                    orthologSets[i] = orthologSets[j];
                    orthologSets[j] = atmp;
                }
            }
        }
        
        for(int i = 0; i< orthologSets.length; i++){
            orthologSets[i].order();
          //  orthologSets[i].print();
        }
        // rename all the genes
        for(int i = 0; i< orthologSets.length; i++){
            for(int j = 0; j< orthologSets[i].genes.length; j++){
                orthologSets[i].genes[j].newName = new Integer(i+1).toString();
            }
            for(int j = 0; j< orthologSets[i].genePairs.length; j++){
                orthologSets[i].genePairs[j].gene1.newName = new Integer(i+1).toString();
                orthologSets[i].genePairs[j].gene2.newName = new Integer(i+1).toString();
            }
        }
        
        

        
    }

    public GenePair[] removeEdge(Homolog ahomolog, int lowSimilarity){
        GenePair[] tmp = new GenePair[ahomolog.genePairs.length];
        int index = 0;
        for(int i = 0; i< ahomolog.genePairs.length; i++){
            if(ahomolog.genePairs[i].weight>lowSimilarity){
                tmp[index] = new GenePair(ahomolog.genePairs[i]);
                index++;
            }
        }
        GenePair[] result = new GenePair[index];
        for(int i = 0; i< result.length; i++){
            result[i] = tmp[i];
        }
        return result;
    }
    public Homolog[] checkGoodOrNot(Homolog ah, int[] allGenomeIndex){
       //ah.print();
		Homolog[] result = new Homolog[1];
		//GeneInfo[] allGenes = ah.getAllGeneList(ah);
		int[][] GenePairMatrix = getGenePairMatrix(ah.genes, ah);
		int numberOfEdges = ah.genePairs.length;
		// public int[] checkGoodMatric(int[][] matrix, GeneInfo[] allGenes, int genomeNumber, int numberOfEdges){
		int[] units= checkGoodMatric(GenePairMatrix, ah.genes, allGenomeIndex,numberOfEdges);
		if(units[2]==0){
			result[0] = ah;
		}else{
			result=new Homolog[0];
		}
		return result;
	}
    
    
    public Homolog[] minimumDeletion(Homolog ah, int[] allGenomeIndex){
        GeneInfo[] allGenes = new GeneInfo[ah.genes.length];
        for(int i = 0; i< allGenes.length; i++){allGenes[i] = new GeneInfo(ah.genes[i]);}
		int[][] GenePairMatrix = getGenePairMatrix(allGenes, ah);
		int[][] GenePairWeightMatrix = getGenePairWeightMatrix(allGenes,ah);
        
		int numberOfEdges = ah.genePairs.length;
        
        GenePairMatrix = matrixAfterDropEdges(GenePairMatrix, GenePairWeightMatrix, allGenes,allGenomeIndex,numberOfEdges);
        
        GenePair[] gps = new GenePair[ah.genePairs.length];
        int index = 0;
        for(int i = 0; i< GenePairMatrix.length; i++){
            for(int j = i+1; j < GenePairMatrix.length; j++){
                if(GenePairMatrix[i][j]>0){
                    gps[index] = new GenePair(allGenes[i],allGenes[j], GenePairWeightMatrix[i][j]);
                    index++;
                }
            }
        }
        GenePair[] newEdges = new GenePair[index];
        for(int i = 0; i< newEdges.length; i++){
            newEdges[i] = gps[i];
        }
        
        Homolog[] result = ah.getAllHomologs1(newEdges,allGenomeIndex.length);
        return result;
	}
    


    public int[][] matrixAfterDropEdges(int[][] GenePairMatrix,int[][] weightMatrix, GeneInfo[] allGenes, int[] allGenomeIndex, int numberOfEdges){
		int[] goodMatrix = checkGoodMatric(GenePairMatrix,allGenes, allGenomeIndex,numberOfEdges);
		Random r = new Random(0);
		while(goodMatrix[2]!=0){
            //	System.out.println("*******bad one****************");
			int besti = -1;
			int bestj = -1;
			int besti2 = -1;
			int bestj2 = -1;
			int biggestGoodChange = -1000;
			int biggestOneChange = 0;
			int smallestWeight = 100;
			for(int i = 0; i< GenePairMatrix.length-1; i++){
				for(int j = i+1; j< GenePairMatrix.length; j++){
					if(GenePairMatrix[i][j]==1){
						int[][] matrixDropAnEdge = new int[GenePairMatrix.length][GenePairMatrix.length];
						for(int k = 0; k<matrixDropAnEdge.length; k++){
							for(int l = 0; l<matrixDropAnEdge.length; l++){
								matrixDropAnEdge[k][l] = GenePairMatrix[k][l];
							}
						}
						matrixDropAnEdge[i][j] = 0;
						matrixDropAnEdge[j][i] = 0;
						int[] goodMatrix2 = checkGoodMatric(matrixDropAnEdge,allGenes,allGenomeIndex,numberOfEdges);
						
						if(goodMatrix2[0]-goodMatrix[0]==biggestGoodChange && goodMatrix[3]-goodMatrix2[3]==biggestOneChange && weightMatrix[i][j]==smallestWeight){
							int replace = r.nextInt(2);
							//	System.out.println("Here" );
							if(replace==1){
                                biggestGoodChange = goodMatrix2[0]-goodMatrix[0];
                                biggestOneChange = goodMatrix[3]-goodMatrix2[3];
                                smallestWeight = weightMatrix[i][j];
                                besti = i;
								bestj = j;}
						}
						
						
						if(goodMatrix2[0]-goodMatrix[0]==biggestGoodChange && goodMatrix[3]-goodMatrix2[3]==biggestOneChange && weightMatrix[i][j]<smallestWeight){
							//	System.out.println("Here" );
							biggestGoodChange = goodMatrix2[0]-goodMatrix[0];
							biggestOneChange = goodMatrix[3]-goodMatrix2[3];
							smallestWeight = weightMatrix[i][j];
							besti = i;
							bestj = j;
						}
						if(goodMatrix2[0]-goodMatrix[0]>biggestGoodChange  || goodMatrix2[0]-goodMatrix[0]==biggestGoodChange && goodMatrix[3]-goodMatrix2[3]>biggestOneChange){
							//	System.out.println("Here" );
							biggestGoodChange = goodMatrix2[0]-goodMatrix[0];
							biggestOneChange = goodMatrix[3]-goodMatrix2[3];
							smallestWeight = weightMatrix[i][j];
							besti = i;
							bestj = j;
						}
					}
				}
			}
			//	System.out.println("============== best i amd best j : "+ besti+ "   "+ bestj);
			if(besti==-1){
				System.out.println("!!!!!!!!!!");
				besti=besti2;
				bestj = bestj2;
			}
			//	System.out.print("d");
			//System.out.print("drop the line "+besti+"  ,  "+ bestj);
			GenePairMatrix[besti][bestj] = 0;
			GenePairMatrix[bestj][besti] = 0;
			goodMatrix = checkGoodMatric(GenePairMatrix,allGenes,allGenomeIndex,numberOfEdges);
			//System.out.println(" :"+ goodMatrix[0]+"-" +goodMatrix[1]+"-"+goodMatrix[2]);
		}
		// add back some good edges
		// order all edges by the weight
		int[][] orderedIndex = getOrderedIndex(weightMatrix);
		for(int i = 0; i< orderedIndex.length; i++){
			if(orderedIndex[i][2]!=0){
                int indexi = orderedIndex[i][0];
                int indexj = orderedIndex[i][1];
                GenePairMatrix[indexi][indexj] = 1;
                GenePairMatrix[indexj][indexi] = 1;
                int[] goodMatrix2 = checkGoodMatric(GenePairMatrix,allGenes,allGenomeIndex,numberOfEdges);
                if(goodMatrix2[2]!=0){
                    GenePairMatrix[indexi][indexj] = 0;
                    GenePairMatrix[indexj][indexi] = 0;
                }
			}
		}
		return GenePairMatrix;
	}

    
    public int[] checkGoodMatric(int[][] matrix, GeneInfo[] allGenes, int[] allGenomeIndex, int numberOfEdges){
        //warshall method
		int[][] finalMatrix = new int[matrix.length][matrix.length];
		for(int i = 0; i< matrix.length; i++){
			for(int j  =0; j < matrix.length; j++){
				finalMatrix[i][j] = matrix[i][j];
			}
		}
		for(int k = 0; k< matrix.length; k++){
			for(int i = 0; i< matrix.length; i++){
				if(finalMatrix[i][k]==1){
					for(int j = 0; j< matrix.length; j++){
						if(finalMatrix[i][j] == 0 && finalMatrix[k][j] == 0){
							finalMatrix[i][j] = 0;
						}else{
							finalMatrix[i][j] = 1;
						}
					}
				}
			}
		}
    
		int badUnit = 0;
		int emptyUnit = 0;
		int goodUnit = 0;
		int ones = 0;
		
		for(int i = 0; i< finalMatrix.length; i++){
			for(int j = 0; j< finalMatrix[0].length; j++){
				if(finalMatrix[i][j] == 1){
					ones++;
				}
			}
		}
		int goodOnes = 0;
		for(int i = 0 ; i< matrix.length; i++){//for each gene
			for(int j = 0; j < allGenomeIndex.length; j++){  //for each genome
                
				int positiveNumber = 0;
				int totalNumber = 0;
				int genekindex = -1;
				for(int k = 0; k< matrix.length; k++){
					if(allGenes[k].genomeIndex==allGenomeIndex[j]){
						genekindex = k;
						totalNumber++;
						if(finalMatrix[i][k]>0){
							positiveNumber++;
						}
					}
				}
				if(genekindex!=-1){
                    boolean findPosition = false;
                    int diploidNumber1= allGenes[i].ploidyNumber;
                    int diploidNumber2 = allGenes[genekindex].ploidyNumber;  // the diploid number of genome j
                    
                    if(findPosition == false && positiveNumber<=diploidNumber2 && positiveNumber>0){
                        goodUnit++;
                        findPosition=true;
                        goodOnes = goodOnes+positiveNumber;
                    }
                    
                    
                    if(findPosition== false && totalNumber>0 && positiveNumber==0 || genekindex==-1){
                        emptyUnit++;
                        findPosition=true;
                    }
                    
                    if(findPosition== false ){
                        badUnit++;
                        findPosition=true;
                    }
				}				
				
			}
		}
		
		int[] result = new int[4];
		result[0] = goodUnit;
		result[1]= emptyUnit;
		result[2] = badUnit;
		result[3] = ones-goodOnes;
		
		return result;
	}
    
    public int[][] getOrderedIndex(int[][] w){
		int[][] indexAndWeight = new int[w.length*w.length][3];
		int index = 0;
		for(int i = 0; i< w.length-1; i++){
			for(int j = i+1; j<w.length; j++){
				if(w[i][j]!=0){
					indexAndWeight[index][0] = i;
					indexAndWeight[index][1] = j;
					indexAndWeight[index][2] = w[i][j];
					index++;
				}
			}
		}
		for(int i = 0; i< indexAndWeight.length-1; i++){
			for(int j = i+1; j<indexAndWeight.length; j++){
				if(indexAndWeight[i][2]<indexAndWeight[j][2]){
					int temp1 = indexAndWeight[i][0];
					int temp2 = indexAndWeight[i][1];
					int temp3 = indexAndWeight[i][2];
					indexAndWeight[i][0] = indexAndWeight[j][0];
					indexAndWeight[i][1] = indexAndWeight[j][1];
					indexAndWeight[i][2] = indexAndWeight[j][2];
					
					indexAndWeight[j][0] = temp1;
					indexAndWeight[j][1] = temp2;
					indexAndWeight[j][2] = temp3;
				}
                
			}
		}
		return indexAndWeight;
	}
    

    
	
    
    public int findTheEntry(GeneInfo[] gl, String ag){
		for(int i = 0; i<gl.length; i++){
			if(gl[i].name.equals(ag)){
				return i;
			}
		}
		System.out.println("something wrong in findTheEntry");
		return -1;
	}
	
	public void printAMatrix(int[][] m){
		for(int i = 0; i< m.length; i++){
			for(int j = 0; j< m.length; j++){
				System.out.print(m[i][j]+"  ");
			}
			System.out.println();
		}
	}
	public int[][] getGenePairMatrix(GeneInfo[] genes, Homolog ah){
		//System.out.println("----"+genes.length);
		int[][] matrix = new int[genes.length][genes.length];
		for(int i = 0; i< ah.genePairs.length; i++){
			if(ah.genePairs[i]!=null){
				String g1 = ah.genePairs[i].gene1.name;
				String g2 = ah.genePairs[i].gene2.name;
				int entry1 = findTheEntry(genes, g1);
				int entry2 = findTheEntry(genes, g2);
				matrix[entry1][entry2] = 1;
				matrix[entry2][entry1] = 1;
			}
		}
		for(int i = 0; i< matrix.length; i++){
			matrix[i][i] = 1;
		}
		return matrix;
	}
	
	
	
	public int[][] getGenePairWeightMatrix(GeneInfo[] genes, Homolog ah){
		//System.out.println("----"+genes.length);
		int[][] matrix = new int[genes.length][genes.length];
		for(int i = 0; i< ah.genePairs.length; i++){
			if(ah.genePairs[i]!=null){
				String g1 = ah.genePairs[i].gene1.name;
				String g2 = ah.genePairs[i].gene2.name;
				int entry1 = findTheEntry(genes, g1);
				int entry2 = findTheEntry(genes, g2);
				matrix[entry1][entry2] = ah.genePairs[i].weight;
				matrix[entry2][entry1] =  ah.genePairs[i].weight;
			}
		}
		for(int i = 0; i< matrix.length; i++){
			matrix[i][i] = 100;
		}
		//	System.out.println("a matrix ");
		//	printAMatrix(matrix);
		return matrix;
	}
	
}
