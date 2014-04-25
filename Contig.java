import java.util.*;
import java.io.*;
public class Contig{
    GenomeInString[] leaveGenomes;
	int contigIndex;
	int totalGeneNumber;
	int[] geneNumber;
	String[] geneName;
    PositionInAGenome[] positionInAGenome;
	int numberOfGenome;
    int[] ploidyNumber;
    int[] allGenomeIndex;
	GeneGroup[] geneGroups;
	int ggIndex;
    GeneGroup[] geneGroupsNoSubgenomeInfo;
	int ggIndexNosi;
    
	
	//public Contig(int numberOfGenome1, int[] allGenomeIndex, int[] ploidyNumber){
//}
	
	
	public Contig(int ci, String agl, int[] ploidyNumber1, GenomeInString[] allGenomes, int[] genomeIndex){
        ploidyNumber = new int[ploidyNumber1.length];
        ploidyNumber = ploidyNumber1;
        
        allGenomeIndex = new int[genomeIndex.length];
        allGenomeIndex = genomeIndex;
        
        int maximumGeneGroup = 0;
        for(int i = 0; i< ploidyNumber.length; i++){
            maximumGeneGroup = maximumGeneGroup+ploidyNumber[i];
        }
		ggIndex = 0;
		geneGroups = new GeneGroup[maximumGeneGroup];
        
        ggIndexNosi = 0;
        geneGroupsNoSubgenomeInfo = new GeneGroup[maximumGeneGroup];
		
		contigIndex = ci; 
		String[] gl = splitBS(agl);
		totalGeneNumber = gl.length;
		geneName = new String[gl.length];
		for(int i = 0; i< gl.length; i++){
			geneName[i] = gl[i];
			if(geneName[i].substring(0,1).equals("-")){
				geneName[i]=geneName[i].substring(1);
			}
		}
       // System.out.println(totalGeneNumber+"\tto get a contig\t"+agl);
        leaveGenomes = new GenomeInString[allGenomes.length];
        leaveGenomes = allGenomes;
        
        positionInAGenome = new PositionInAGenome[ploidyNumber.length];
        for(int i = 0; i< positionInAGenome.length; i++){
            positionInAGenome[i] = new PositionInAGenome(ploidyNumber[i],totalGeneNumber);
        }
        
        for(int i = 0; i< geneName.length; i++){
            for(int j  =0; j< ploidyNumber.length; j++){
                int[][] position1 = findGenePosition(geneName[i],leaveGenomes[j],ploidyNumber[j]);
                positionInAGenome[j].addPositionForAGene(i, position1);
            }
        }
		
		geneNumber = new int[gl.length];
        
      /*  for(int i = 0; i< ploidyNumber.length; i++){
            positionInAGenome[i].print();
        }*/
		
	}
    
    public Contig(int ci, GeneGroup agg){
        contigIndex = ci;
        geneGroupsNoSubgenomeInfo = new GeneGroup[1];
        geneGroupsNoSubgenomeInfo[0] = new GeneGroup(agg);
        ggIndexNosi=1;
    }
	
	
	public int[][] findGenePosition(String ag, GenomeInString agenome, int d){
		int[][] result = new int[d][2];for(int i = 0; i< d; i++){for(int j = 0; j<2; j++){result[i][j]=-1;}}
		int index = 0;
		//System.out.println(ag);
		for(int i = 0; i< agenome.chrs.length; i++){
			String[] genes = splitBS(agenome.chrs[i]);
			for(int j = 0; j< genes.length; j++){
				String aghere = genes[j];
				if(aghere.substring(0,1).equals("-")){
					aghere = aghere.substring(1);
				}
				if(ag.equals(aghere)){
					result[index][0] = i;
					result[index][1] = j;
					index++;
				}
			}
		}
		return result;
	}
    
    
	
	public void printAContig(){
		System.out.println("======================="+contigIndex+"\t"+totalGeneNumber);
		System.out.print("geneName\t");
		for(int i = 0; i< geneName.length; i++){
			System.out.print(geneName[i]+"\t");
		}
		System.out.println();
        for(int i = 0; i< ploidyNumber.length; i++){
            for(int j = 0; j< ploidyNumber[i]; j++){
                System.out.print("genome"+allGenomeIndex[i]+"-copyChr"+j+"\t");
                for(int k = 0; k< geneName.length; k++){
                    System.out.print(positionInAGenome[i].chrs[j][k]+"\t");
                }
                System.out.println();
                System.out.print("genome"+allGenomeIndex[i]+"-copyPosition"+j+"\t");
                for(int k = 0; k< geneName.length; k++){
                    System.out.print(positionInAGenome[i].positions[j][k]+"\t");
                    
                }
                 System.out.println();
            }
        }
        for(int i = 0; i<ggIndex; i++){
			geneGroups[i].printAGeneGroup();
		}
        for(int i = 0; i<ggIndexNosi; i++){
			geneGroupsNoSubgenomeInfo[i].printAGeneGroup();
		}
        
	}
    
    public int[] getPositionList(int c, int[][] chrs, int[][] pos){
		int[] tmp = new int[pos.length*pos[0].length];
		int index = 0;
		for(int i = 0; i< chrs.length; i++){
			for(int j = 0; j< chrs[i].length; j++){
				if(chrs[i][j]==c){
					tmp[index] = pos[i][j];
					index++;
				}
			}
		}
		int[] result = new int[index];
		for(int i = 0; i< index; i++){
			result[i]= tmp[i];
		}
		/*for(int i = 0; i< result.length-1; i++){
         for(int j = i+1; j<result.length; j++){
         if(result[i]>result[j]){
         int atmp = result[i];
         result[i] = result[j];
         result[j] = atmp;
         }
         }
        }*/
		return result;
	}
    
    
    public int[] getMainChr(int[][] chr, int diploidNumber){
		int[] result = new int[diploidNumber];for(int i = 0; i< result.length; i++){result[i]=-1;}
        int biggestChrNumber = 0;
        for(int i = 0; i< chr.length; i++){
            for(int j = 0; j< chr[0].length; j++){
                if(chr[i][j]> biggestChrNumber){
                    biggestChrNumber = chr[i][j];
                }
            }
        }
		int[] count = new int[biggestChrNumber+5];
		for(int i = 0; i< chr.length; i++){
			for(int j = 0; j< chr[i].length; j++){
                if(chr[i][j]!=-1){
                    count[chr[i][j]]++;
                }}
		}
		for(int i = 0; i< result.length; i++){
			int biggest = 0;
			int chrHere = -1;
			for(int j = 0; j< count.length; j++){
				if(count[j]>biggest){
					biggest = count[j];
					chrHere = j;
				}
			}
			if(chrHere!=-1){
                result[i] = chrHere;
				count[chrHere] = 0;}
		}
		return result;
	}

				
	public void getGeneGroupInAGenome(int minimumSize, int gi, int[][] biRange){
        int ploidy = ploidyNumber[gi];
      //  positionInAGenome[gi].print();
        int[][] chr = positionInAGenome[gi].chrs;
        int[][] pos = positionInAGenome[gi].positions;
		
		int[] mainChr= getMainChr(chr,ploidy);
		//System.out.println("mainChr :"+ mainChr[0]+"  "+mainChr[1]+"  "+ mainChr[2]);
		for(int i = 0; i< mainChr.length; i++){
			if(mainChr[i]!=-1){
			//	System.out.println("when mainChr = "+ mainChr[i]);
			int[] posList = getPositionList(mainChr[i], chr, pos);
				//System.out.print("posList ");for(int j = 0; j< posList.length; j++){System.out.print(posList[j]+"  ");}System.out.println();
				
			int[] goodOne = checkGoodOrNot(mainChr[i],posList,gi, minimumSize, biRange);
				
			if(goodOne[0] !=-1){
          //  if(goodOne.length!=0){
			//	System.out.print("posList ");for(int j = 0; j< posList.length; j++){System.out.print(posList[j]+"  ");}System.out.println();
				geneGroups[ggIndex]= retriveAGeneGroup(gi, goodOne, mainChr[i], chr, pos);
				ggIndex++;
			}}
		}
	}
    
    public GeneGroup retriveAGeneGroup(int gi, int[] values, int c, int[][] chr, int[][] pos ){
		GeneGroup result = new GeneGroup();
		result.genomeIndex = gi;
		result.geneNumber = values.length-2;
		result.geneName = new String[result.geneNumber];
        //System.out.println("geneNumber ")
      //  System.out.println("chrLenght "+ chr.length+"\t"+ chr[0].length);
		for(int i = 2; i<values.length; i++){
			int poshere = values[i];
			for(int j = 0; j< chr.length; j++){
				for(int k = 0; k< chr[j].length; k++){
					if(chr[j][k]==c && pos[j][k]==poshere){
                      //  System.out.println("here "+ i+"\t"+j);
						result.geneName[i-2] = geneName[k];
					}
				}
			}
		}
		result.chr = c;
		result.pos = new int[values.length-2];
		for(int i = 0; i< values.length-2; i++){
			result.pos[i] =  values[i+2];
		}
		double sum = 0;
		for(int i = 0; i< result.pos.length; i++){
			sum = (double)(sum+result.pos[i]);
		}
		result.averagePos = (double)sum/result.pos.length;
		result.colorCode = values[0];
		result.blockIndex = values[1];
		return result;
	}

    
    public void getGeneGroupNoSubgenomeInfoInAGenome(int minimumSize, int gi, int[][] biRange){
        int ploidy = ploidyNumber[gi];
        int[][] chr = positionInAGenome[gi].chrs;
        int[][] pos = positionInAGenome[gi].positions;
		int[] mainChr= getMainChr(chr,ploidy);
		//System.out.println("mainChr :"+ mainChr[0]+"  "+mainChr[1]+"  "+ mainChr[2]);
		for(int i = 0; i< mainChr.length; i++){
			if(mainChr[i]!=-1){
				//System.out.println("when mainChr = "+ mainChr[i]);
                int[] posList = getPositionList(mainChr[i], chr, pos);
				//System.out.print("posList ");for(int j = 0; j< posList.length; j++){System.out.print(posList[j]+"  ");}System.out.println();
                if(posList.length >=minimumSize){
                    GeneGroup result = new GeneGroup();
                    result.genomeIndex = gi;
                    result.geneNumber = posList.length;
                    result.geneName = new String[result.geneNumber];
                    
                    for(int p = 0; p<posList.length; p++){
                        int poshere = posList[p];
                        for(int j = 0; j< chr.length; j++){
                            for(int k = 0; k< chr[j].length; k++){
                                if(chr[j][k]==mainChr[i] && pos[j][k]==poshere){
                                    result.geneName[p] = geneName[k];
                                }
                            }
                        }
                    }
                    result.chr = mainChr[i];
                    result.pos = new int[posList.length];
                    for(int p = 0; p< posList.length; p++){
                        result.pos[p] =  posList[p];
                    }
                    double sum = 0;
                    for(int p = 0; p< result.pos.length; p++){
                        sum = (double)(sum+result.pos[p]);
                    }
                    result.averagePos = (double)sum/result.pos.length;
                    result.colorCode = -1;
                    result.blockIndex = -1;
                    
                    geneGroupsNoSubgenomeInfo[ggIndexNosi] = result;
                    ggIndexNosi++;
                }
            }
		}
	}
    
    
    public int[] findBIIndex(int[][] biRange, int c, int[] posList, int order){
		int[] result = new int[3]; result[0]=-1; result[1] = -1; result[2] = -1;
		for(int i = 0; i< posList.length; i++){
			int pos = posList[i];
			if(order == 1){pos = posList[posList.length-1-i];}
			for(int j = 0; j< biRange.length; j++){
				if(biRange[j][2]==c && biRange[j][3]<=pos && biRange[j][4]>=pos){
					result[0] = biRange[j][0];
					result[1] = biRange[j][1];
					result[2] = i;
					if(order==1){result[2] = posList.length-1-i;}
					return result;
				}
			}
		}
		return result;
	}
	
    
	
   //***************
    //**************
    // if biRange.length ==0;
	public int[] checkGoodOrNot(int mainChr, int[] posList, int genomeIndex, int minimumSize,  int[][] biRange){
		int[] result = new int[posList.length]; for(int i = 0; i< result.length; i++){result[i] = -1;}
		if(posList.length<minimumSize){result[0] = -1; return result;
          //  result = new int[0];
           // return result;
        }
       
        
		int[] firstBiIndex = findBIIndex(biRange, mainChr,posList, 0);
	//	System.out.println("firstBiIndex "+ firstBiIndex[0]+"  "+ firstBiIndex[1] +"  "+ firstBiIndex[2]);
		int[] lastBiIndex = findBIIndex(biRange, mainChr,posList, 1);
	//	System.out.println("lastBiIndex "+ lastBiIndex[0]+"  "+ lastBiIndex[1] +"  "+ lastBiIndex[2]);

		if(firstBiIndex[1]!=lastBiIndex[1]) {
			int lastFirstBI = findTheOtherEnd(mainChr, posList, firstBiIndex[1], 0, biRange);
			int firstLastBI = findTheOtherEnd(mainChr, posList, lastBiIndex[1],1,biRange);
		//	System.out.println("lastFirstBI  and firstLastBI = "+ lastFirstBI+"  "+firstLastBI);
			int l1 = lastFirstBI-firstBiIndex[2]+1;
			int l2 = lastBiIndex[2]-firstLastBI +1;
		//	System.out.println("l1  and l2 = "+ l1+"  "+l2);
			if(l1>l2 && l1>= minimumSize){
				result = new int[l1+2];
				result[0] = firstBiIndex[0];
				result[1] = firstBiIndex[1];
				for(int i = 2; i< result.length; i++){
					result[i] = posList[i+firstBiIndex[2]-2];
				}
			}
			if(l2>l1 && l2>=minimumSize){
				result = new int[l2+2];
				result[0] = lastBiIndex[0];
				result[1] = lastBiIndex[1];
				
				for(int i = 2; i< result.length; i++){
					result[i] = posList[i+firstLastBI-2];
				}
			}
		}
		if(firstBiIndex[1]==lastBiIndex[1] && firstBiIndex[0]!=-1){
			int l = (lastBiIndex[2]-firstBiIndex[2]+1);
			if(l>minimumSize){
				result = new int[lastBiIndex[2]-firstBiIndex[2]+1+2];
				result[0] = firstBiIndex[0];
				result[1] = firstBiIndex[1];
			//System.out.println("l1 "+ (lastBiIndex[2]-firstBiIndex[2]+1));
				for(int i = 2; i< result.length; i++){
					result[i] = posList[i+firstBiIndex[2]-2];
				}
			}
			
		}
        
		//System.out.println("values " );
		//for(int i = 0; i< result.length; i++){System.out.print(result[i]+"  ");}System.out.println();
		return result;
	}
	
	public int findTheOtherEnd(int c, int[] pl, int bi, int order, int[][] biRange){
		for(int i = 0; i< pl.length; i++){
			int poshere = pl[pl.length-1-i];
			if(order == 1){ poshere = pl[i];}
			for(int j = 0; j< biRange.length; j++){
				if(biRange[j][2]==c && biRange[j][3]<=poshere && biRange[j][4]>=poshere){
					if(biRange[j][1] == bi){
						if(order ==0){return pl.length-1-i;}
						if(order == 1){return i;}
					}
				}
			}
		}
		return -1;
	}
    
    public String[] splitBS(String s){
        String[] tmp = s.trim().split(" ");
        int index = 0;
        for(int i = 0; i < tmp.length; i++){
            if(tmp[i].trim().equals("") ==false){
                index++;
            }
        }
        String[] result = new String[index];
        index = 0;
        for(int j = 0; j <tmp.length; j++){
            if(tmp[j].trim().equals("") == false){
                result[index] = tmp[j].trim();
                index++;
            }
        }
        return result;
    }
    
    
	
	   
    public boolean inColorCode(Contig ac, int cc){
        /*if(ac.contigIndex==476){
            ac.printAContig();
            for(int i = 0; i< ac.ggIndex; i++){
                System.out.println("cchere "+ac.geneGroups[i].colorCode);
            }
        }*/
		if(ac.ggIndex==0){return false;}
		for(int i = 0; i< ac.ggIndex; i++){
			if(ac.geneGroups[i].colorCode!=cc){
				return false;
			}
		}
		return true;
	}

    
    // gene content
	public GenomeInString getGenomeGeneContent(Contig[] contigs){
		String[] genome = new String[1];genome[0]= "";
		for(int i = 0; i< contigs.length; i++){
			if(contigs[i]!=null){
				genome[0]= genome[0] +"  contig"+new Integer(contigs[i].contigIndex).toString();
			}
		}
		GenomeInString ag = new GenomeInString(genome);
		return ag;
	}
    
	
	public GenomeInString[] getGenomesInContig(Contig[] allContig, int colorCode, SubgenomeRanges[] rangesForGenomes){
        Contig[] tmp  = new Contig[allContig.length];
        int genomeNumberhere =1;
        for(int i = 0;i< ploidyNumber.length; i++){
            if(rangesForGenomes[i].ranges.length==0){
                genomeNumberhere = genomeNumberhere+1;
            }else{
                genomeNumberhere= genomeNumberhere+ploidyNumber[i];
            }
        }
        if(colorCode!=-1){
            int tmpIndex = 0;
            for(int i = 0; i< allContig.length; i++){
                if(allContig[i]!=null && inColorCode(allContig[i],colorCode)){
                    tmp[tmpIndex] = allContig[i];
                    tmpIndex++;
                }
            }
        }else{
            int tmpIndex = 0;
            for(int i = 0; i< allContig.length; i++){
                if(allContig[i]!=null && allContig[i].ggIndexNosi!=0){
                    tmp[tmpIndex] = allContig[i];
                    tmpIndex++;
                }
            }
           // tmp = allContig;
        }
        
		//System.out.println("contigs in ancestor of colorCode "+colorCode);
		GenomeInString[] allGenomesInContig = new GenomeInString[genomeNumberhere];
		allGenomesInContig[0] = getGenomeGeneContent(tmp);  // gene content
        int index = 1;
        
        for(int i = 0; i< leaveGenomes.length; i++){
            if(rangesForGenomes[i].ranges.length==0){
                Contig[] contigInAGenome = orderContigNoSubgenomeInfoForAGenome(tmp, i);
                allGenomesInContig[index] = getGenomeNoSubgenomeInfo(contigInAGenome,i);
                index++;
            }else{
            for(int j = 0; j< ploidyNumber[i]; j++){
                Contig[] contigInAGenome = orderContigByGeneGroupForAGenome(tmp, colorCode, j, i);
                allGenomesInContig[index] = getGenome(contigInAGenome, colorCode,j,i);
                index++;
            }}
        }
        
		/*for(int i = 0; i< allGenomesInContig.length; i++){
            System.out.println(i+"  which genome");
			int gn = countGeneNumber(allGenomesInContig[i].chrs);
			System.out.println("genome "+i+ "   geneNumber= "+gn);
			for(int j = 0; j< allGenomesInContig[i].chrs.length; j++){
				System.out.println("chr "+j+"\n"+allGenomesInContig[i].chrs[j]);
			}
		}*/
		
		return allGenomesInContig;
		
		
	}
    
    public int getGeneGroupIndex(Contig ac, int cc, int bi, int gi){
		for(int i = 0; i< ac.ggIndex; i++){
			if(ac.geneGroups[i].colorCode == cc && ac.geneGroups[i].blockIndex == bi && ac.geneGroups[i].genomeIndex == gi){
				return i;
			}
		}
		return -1;
	}
	

    public String checkSignGeneGroup(int[] p){
		int inc = 0;
		int dec = 0;
		int startIndex = 0;
		int pre = -1;
		pre = p[0];
		for(int i = 1; i< p.length; i++){
            if(p[i]>pre){inc++;}
            if(p[i]<pre){dec++;}
            pre = p[i];
		}
		//	System.out.println("incre "+ inc);
		//	System.out.println("dec "+ dec);
		
		if(inc>=dec){
			return "";
		}
		return "-";
	}

    public GenomeInString getGenomeNoSubgenomeInfo(Contig[] contigs, int gi){
		int biggestChr = 0;
        for(int i = 0; i< leaveGenomes.length; i++){
            if(leaveGenomes[i].chrs.length> biggestChr){biggestChr = leaveGenomes[i].chrs.length;}
        }
		String[] genome = new String[biggestChr];for(int i = 0; i< genome.length; i++){ genome[i] = "";}
		for(int i = 0; i< contigs.length; i++){
			if(contigs[i]!=null){
                String sign = checkSignGeneGroup(contigs[i].geneGroupsNoSubgenomeInfo[0].pos);
                int chr = contigs[i].geneGroupsNoSubgenomeInfo[0].chr;
                genome[chr]= genome[chr] +"  "+sign+"contig"+new Integer(contigs[i].contigIndex).toString();
				
			}
        }
    
		// remove empty chr
		int emptyChr = 0;
		for(int i = 0; i< genome.length; i++){
			if(genome[i].equals("")){emptyChr++;}
		}
		String[] result = new String[biggestChr-emptyChr];
		int index = 0;
		for(int i = 0; i< genome.length; i++){
			if(genome[i].equals("")==false){
				result[index] = genome[i];
				index++;
			}
		}
		GenomeInString ag = new GenomeInString(result);
		return ag;
    }

    
    public GenomeInString getGenome(Contig[] contigs, int cc, int bii, int gi){
		int biggestChr = 0; for(int i = 0; i< leaveGenomes.length; i++){
            if(leaveGenomes[i].chrs.length> biggestChr){biggestChr = leaveGenomes[i].chrs.length;}
        }
		String[] genome = new String[biggestChr];for(int i = 0; i< genome.length; i++){ genome[i] = "";}
		for(int i = 0; i< contigs.length; i++){
			if(contigs[i]!=null){
				int indexi1 = getGeneGroupIndex(contigs[i], cc,bii, gi);
				if(indexi1!=-1){
					String sign = checkSignGeneGroup(contigs[i].geneGroups[indexi1].pos);
					int chr = contigs[i].geneGroups[indexi1].chr;
					genome[chr]= genome[chr] +"  "+sign+"contig"+new Integer(contigs[i].contigIndex).toString();
				}
			}
		}
		// remove empty chr
		int emptyChr = 0;
		for(int i = 0; i< genome.length; i++){
			if(genome[i].equals("")){emptyChr++;}
		}
		String[] result = new String[biggestChr-emptyChr];
		int index = 0;
		for(int i = 0; i< genome.length; i++){
			if(genome[i].equals("")==false){
				result[index] = genome[i];
				index++;
			}
		}
		GenomeInString ag = new GenomeInString(result);
		return ag;
    }

    public int countGeneNumber(String[] g1){
		int gn2 = 0;
		for(int i = 0; i< g1.length; i++){
			String[] genes = splitBS(g1[i]);
			gn2 = gn2+genes.length;
		}
		//System.out.println("gene number ="+ gn2);
		return gn2;
	}
  //  GeneGroup[] geneGroupsNoSubgenomeInfo;
	//int ggIndexNosi;

    public Contig[] orderContigNoSubgenomeInfoForAGenome(Contig[] contigs, int gi){ // no subgenome infomation
       /* System.out.println("before order contig no subgenomeInf0");
        for(int i = 0; i< contigs.length; i++){
            if(contigs[i]!=null){System.out.println("contig "+contigs[i].contigIndex);}
        }*/
        int ploidyhere = ploidyNumber[gi];
        Contig[] tmp = new Contig[contigs.length*ploidyhere];
        int tmpIndex = 0;
        for(int i = 0; i< contigs.length; i++){
            if(contigs[i]!=null){
            for(int j = 0; j< contigs[i].geneGroupsNoSubgenomeInfo.length; j++){
                if(contigs[i].geneGroupsNoSubgenomeInfo[j]!=null){
                GeneGroup gghere =contigs[i].geneGroupsNoSubgenomeInfo[j];
                if(gghere!=null && gghere.genomeIndex == gi){
                    tmp[tmpIndex] = new Contig(contigs[i].contigIndex, gghere);
                    tmpIndex++;
                }}
            }}
        }
    //    System.out.println("total contig "+ tmpIndex);
        Contig[] result = new Contig[tmpIndex];
        for(int i = 0; i< result.length; i++){
            result[i] = tmp[i];
        }
        for(int i = 0; i< result.length-1; i++){
            for(int j = i+1; j< result.length; j++){
                GeneGroup gg1 = result[i].geneGroupsNoSubgenomeInfo[0];
                GeneGroup gg2 = result[j].geneGroupsNoSubgenomeInfo[0];
                if(gg1.chr> gg2.chr || gg1.chr==gg2.chr && gg1.averagePos> gg2.averagePos){
                    Contig atmp = result[i];
                    result[i] = result[j];
                    result[j] = atmp;
                }
            }
        }
     /*   for(int i = 0; i< result.length; i++){
            System.out.print(result[i].contigIndex+"  ");
        }
        System.out.println();*/

        return result;
        
    }

    public Contig[] orderContigByGeneGroupForAGenome(Contig[] contigs,  int cc, int bii, int gi){
        
		Contig[] tmp = new Contig[contigs.length];
		int index = 0;
		for(int i = 0; i< contigs.length; i++){
			if(contigs[i]!=null){
				for(int j = 0; j< contigs[i].ggIndex; j++){
					if(contigs[i].geneGroups[j].colorCode == cc && contigs[i].geneGroups[j].blockIndex == bii && contigs[i].geneGroups[j].genomeIndex == gi){
						tmp[index] = contigs[i];
						index++;
						break;
					}
				}
			}
		}
		contigs = new Contig[index];
		for(int i = 0; i< contigs.length; i++){
			contigs[i] = tmp[i];
		}
      /*  System.out.println("======================cc and bi "+ cc +" " + bii+"  genomeIndex "+ gi);
        for(int i = 0; i< contigs.length; i++){
            System.out.print(contigs[i].contigIndex+"  ");
        }
        System.out.println();*/
        
        for(int i = 0; i< contigs.length-1; i++){
            for(int j = i+1; j< contigs.length; j++){
                if(contigs[i]!=null && contigs[j]!=null){
                    int indexi = getGeneGroupIndex(contigs[i], cc,bii, gi);
                    int indexj = getGeneGroupIndex(contigs[j],cc,bii,gi);
                    int chri=-1; int chrj = -1; double mpi = -1;  double mpj = -1;
                    if(indexi!=-1){chri = contigs[i].geneGroups[indexi].chr; mpi = calAverage(contigs[i].geneGroups[indexi].pos);}
                    if(indexj!=-1){chrj = contigs[j].geneGroups[indexj].chr; mpj = calAverage(contigs[j].geneGroups[indexj].pos);}
                    if(chri==-1 || (chri>chrj && chrj!=-1)  || (chri!=-1 && chri== chrj && mpi>mpj)){
                        Contig atmp = contigs[i];
                        contigs[i] = contigs[j];
                        contigs[j] = atmp;
                    }
                }
            }
        }
         
       /*  System.out.println("======================cc and bi   after order "+ cc +" " + bii+"  genomeIndex "+ gi);
         for(int i = 0; i< contigs.length; i++){
             System.out.print(contigs[i].contigIndex+"  ");
         }
         System.out.println();*/
        return contigs;
         
    }

	
		
	public double calAverage(int[] p){
		double total = 0;
		for(int i = 0; i< p.length; i++){
			total = total+p[i];
		}
		return total/p.length;
	}
	
	

	
		
	
}