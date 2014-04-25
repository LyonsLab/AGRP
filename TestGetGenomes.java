
import java.util.*;
import java.io.*;
import java.util.regex.Pattern;


public class TestGetGenomes{
	public static void main( String[] args) throws Exception{
		
		/*Calendar time = Calendar.getInstance();
		System.out.println("start " + time.getTime());
*/
        TestGetGenomes atest = new TestGetGenomes();
        if (args.length < 8) {
            System.out.println(Usage());
            System.exit(1);
        }
        
        String[] genomeIndex = null;
        String[] ploidyArray = null;
        String rangeFile = null;
        File[] data = null;
        String output_directory = ".";
        
        for(int i = 0; i < args.length; i += 2) {
            String argument = args[i].substring(1);
            String option = args[i + 1];
            
            switch(argument) {
                case "d":
                    String dir = option;
                    File f = new File(dir);
                    data = f.listFiles();
                    break;
                case "g":
                    String genomeIndexString = option;
                    genomeIndex = option.split(Pattern.quote(","));
                    break;
                case "p":
                    String ploidyString = option;
                    ploidyArray = ploidyString.split(Pattern.quote(","));
                    break;
                case "s":
                    File sub_genome_file = new File(option);
                    if (!sub_genome_file.exists()) {
                        System.out.println("The sub genome file does not exist.");
                        System.exit(1);
                    }
                    rangeFile = option;
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
        
        if (genomeIndex == null || rangeFile == null || ploidyArray == null ||
            rangeFile == null || data == null) {
            System.out.println("Please specify all the command line arguments");
            System.out.println(Usage());
            System.exit(1);
        }
        
        File fd = new File(output_directory);
        if (!fd.exists()) {
            fd.mkdir();
        }
        
        //int[] genomeInde
        int numberOfGenomes = genomeIndex.length;
        int[] allGenomeIndex = new int[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            allGenomeIndex[i] = Integer.parseInt(genomeIndex[i]);
        }
        
        int[][] ploidy = new int[numberOfGenomes][2];
        for(int i = 0; i< ploidy.length;i++){
            ploidy[i][0] = allGenomeIndex[i];
            ploidy[i][1] =Integer.parseInt(ploidyArray[i]);
            System.out.println(ploidy[i][0]+"\t"+ploidy[i][1]);
        }
        
        
        SubgenomeRanges[] rangesForGenomes = new SubgenomeRanges[numberOfGenomes];
        FileReader fr = new FileReader(rangeFile);
        BufferedReader br = new BufferedReader(fr);
        String aline = "";
        for(int g = 0; g< rangesForGenomes.length; g++){
            aline = br.readLine();
            String[] info = aline.split("\t");
            int gi = Integer.parseInt(info[0]);
            int b = Integer.parseInt(info[1]);
            int p = Integer.parseInt(info[2]);
            System.out.println(gi +"\t"+b+"\t"+p);
            int[][] r = new int[b][5];
            aline = br.readLine();
            for(int i = 0; i< r.length; i++){
                aline = br.readLine();
                info = aline.split("\t");
                for(int j = 0; j< info.length; j++){
                    r[i][j] = Integer.parseInt(info[j]);
                }
            }
            rangesForGenomes[g] =new SubgenomeRanges(gi, p,b, r);
        }
    
     
        int numberOfGenomePairs = data.length;
        GenomePair[]  gps = new GenomePair[numberOfGenomePairs];
        for(int i = 0; i< numberOfGenomePairs; i++){
           // String[] apair = inputFiles[i].split("\t");
            File aFile = data[i];
           // String aFileName =  apair[2];
           // System.out.println(aFileName);
        
            int[] twoGenomeIndex = atest.getGenomeIndex(aFile.getName()); //
            int genomeIndex1 = twoGenomeIndex[0];
            int genomeIndex2 = twoGenomeIndex[1];
            int ploidy1 = atest.getPloidNumber(genomeIndex1, ploidy);
            int ploidy2 = atest.getPloidNumber(genomeIndex2, ploidy);
            System.out.println(genomeIndex1+"\t"+genomeIndex2+"\t"+ploidy1+"\t"+ploidy2);
            
            fr = new FileReader(aFile);
            br = new BufferedReader(fr);
            int totalLines = 0;
            aline = br.readLine();
            while(aline!=null){
                totalLines++;
                aline = br.readLine();
            }
            System.out.println(totalLines);
            
            GeneInfo[] geneList1 = new GeneInfo[totalLines];
            GeneInfo[] geneList2 = new GeneInfo[totalLines];
            int[] weight = new int[totalLines];
            int index = 0;
            fr = new FileReader(aFile);
            br = new BufferedReader(fr);
            aline = br.readLine();
            while(aline!=null){
                if(aline.equals("#")==false && aline.length()>200){
                    String[] infos = aline.split("\t");
                    String[] g1info = atest.splitByTwoLine(infos[1]);
                    String gn1=g1info[3];
                    int po1start =new Integer(g1info[1]).intValue();
                    int po1end = new Integer(g1info[2]).intValue();
                    int chr1 = atest.getChrNumber(g1info[0], genomeIndex1);
                    String sign1="+"; if(g1info[4].equals("-1")){sign1="-";}
                    String[] g2info = atest.splitByTwoLine(infos[5]);
                    String gn2=g2info[3];
                    int po2start = new Integer(g2info[1]).intValue();
                    int po2end =new Integer(g2info[2]).intValue();
                    int chr2 = atest.getChrNumber(g2info[0],genomeIndex2);
                    String sign2="+";if(g2info[4].equals("-1")){sign2="-";}
                    double weighthere = Double.parseDouble(g1info[8]);
                    int wi = (int)(weighthere);
                    if(chr1!=-1 && chr2!=-1 && atest.notInList(geneList1, geneList2, gn1, gn2)){
                        weight[index] = (int)(weighthere);
                       // weightDistribution[wi]++;
                        geneList1[index] = new GeneInfo();
                        geneList1[index].name = gn1;
                        geneList1[index].chrNumber = chr1;
                        geneList1[index].positionStart = po1start;
                        geneList1[index].positionEnd = po1end;
                        geneList1[index].orientation = sign1;
                        geneList1[index].genomeIndex = genomeIndex1;
                        geneList1[index].ploidyNumber = ploidy1;
					
                        geneList2[index] = new GeneInfo();
                        geneList2[index].name = gn2;
                        geneList2[index].chrNumber = chr2;
                        geneList2[index].positionStart = po2start;
                        geneList2[index].positionEnd = po2end;
                        geneList2[index].orientation = sign2;
                        geneList2[index].newName = (new Integer(index+1)).toString();
					
                        geneList2[index].genomeIndex = genomeIndex2;
                        geneList2[index].ploidyNumber = ploidy2;
                        index++;
                       // System.out.print(index+" ");
                    }
                }
                 aline = br.readLine();
              //  System.out.println(aline);
            }
            gps[i] = new GenomePair(geneList1, geneList2,weight, index);
        }
    
        Homolog ahomolog = new Homolog();
		Homolog[] allOriginalHomologs = ahomolog.getAllHomologs(gps, numberOfGenomes);
       
        int[] homologSize1 = new int[10000];
        int[] homologSize2 = new int[10000];
        
        for(int i = 0; i< allOriginalHomologs.length; i++){
            homologSize1[allOriginalHomologs[i].genePairs.length]++;
            if(allOriginalHomologs[i].genePairs.length==38){
              //  allOriginalHomologs[i].print();
            }
        }
        
        for(int i = 0; i< allOriginalHomologs.length; i++){
            homologSize2[allOriginalHomologs[i].genes.length]++;
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
        
        for(int i = 0; i< allOriginalHomologs.length-1; i++){
            for(int j = i+1; j< allOriginalHomologs.length; j++){
                for(int k = 0; k< allOriginalHomologs[i].genes.length; k++){
                    for(int h = 0; h< allOriginalHomologs[j].genes.length; h++){
                        if(allOriginalHomologs[i].genes[k].name.equals(allOriginalHomologs[j].genes[h].name)){
                            System.out.println("more combine "+i+"\t"+j);
                        }
                    }
                }
            }
        }
     
        
        // method from: OMG! Orthologs in multiple genomes - competing graph-theoretical formulations (C. Zheng, K.M. Swenson, E. Lyons, & D. Sankoff),WABI'11 Proceedings of the 11th International Conference on Algorithms in Bioinformatics (T. Przytycka & MF. Sagot ed.) Lecture Notes in Computer Science 6833, 364-375 (2011).
         // Method can be replaced with others for constraining gene families.
		OMGMEC omgmec  = new OMGMEC(allOriginalHomologs, allGenomeIndex);
        omgmec.MECmethod(500,80);
        
        
        // print out orthologSet list in a file
        String outPutFileName = output_directory + "/orthologSets";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        outPutFileName = outPutFileName+".txt";
        FileWriter fstream = new FileWriter(outPutFileName,false);
        BufferedWriter fbw = new BufferedWriter(fstream);

        
        for(int i = 0; i< omgmec.orthologSets.length; i++){
            fbw.write((i+1)+"\t");
            fbw.write(omgmec.orthologSets[i].toStringOrthologSet(allGenomeIndex)+"\n");
        }
        fbw.flush();
        fbw.close();

    
		GenomeInGeneInfo[] genomes = new GenomeInGeneInfo[numberOfGenomes];
        
        outPutFileName = output_directory + "/genomesInString";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        outPutFileName = outPutFileName+".txt";
        fstream = new FileWriter(outPutFileName,false);
        fbw = new BufferedWriter(fstream);
        fbw.write("totalGeneNumber\t"+ omgmec.orthologSets.length+"\n");
        GenomeInString[] genomesInString = new GenomeInString[numberOfGenomes];
		for(int i = 0; i< genomes.length; i++){
            genomes[i] = new GenomeInGeneInfo(omgmec.orthologSets,allGenomeIndex[i]);
			System.out.println("geneNumber in genome\t"+ allGenomeIndex[i] +":\t"+genomes[i].genes.length);
			genomes[i].orderGene();
			genomesInString[i] = genomes[i].getGenomeInString();
            fbw.write("genomeIndex\t"+allGenomeIndex[i]+"\tchrNumber\t"+ genomesInString[i].chrs.length+"\n");
        }
        
        for(int i = 0; i< genomesInString.length; i++){
            fbw.write("genomeIndex\t"+allGenomeIndex[i]+"\n");
            for(int j = 0; j< genomesInString[i].chrs.length; j++){
                fbw.write("chr\t"+j+"\n");
                fbw.write(genomesInString[i].chrs[j]+"\n");
            }
		}
        fbw.flush();
        fbw.close();
        
        //convert ranges in basepairs into geneorders
        
        SubgenomeRanges[] rangesForGenomesInGeneOrder = new SubgenomeRanges[numberOfGenomes];
        for(int i = 0; i< numberOfGenomes; i++){
            rangesForGenomesInGeneOrder[i] = rangesForGenomes[i].getRangesInGeneOrder(genomes[i]);
        }
        
        outPutFileName = output_directory + "/subgenomeRangesInGeneOrder";
        for(int i = 0; i< allGenomeIndex.length; i++){
            outPutFileName = outPutFileName+"_"+new Integer(allGenomeIndex[i]).toString();
        }
        outPutFileName = outPutFileName+".txt";
        fstream = new FileWriter(outPutFileName,false);
        fbw = new BufferedWriter(fstream);
        for(int i = 0; i< rangesForGenomesInGeneOrder.length; i++){
            String rangesForAGenomeInString = rangesForGenomesInGeneOrder[i].toString();
            fbw.write(rangesForAGenomeInString);
        }
        fbw.flush();
        fbw.close();
        
        
	}
	

	public TestGetGenomes(){
	}
    
    public int getPloidNumber(int aGenomeIndex, int[][] ploidyList){
        for(int i = 0; i< ploidyList.length; i++ ){
            if(aGenomeIndex == ploidyList[i][0]){return ploidyList[i][1];}
        }
        return -1;
    }

    public String[] splitByTwoLine(String as){
        String[] tmp = new String[as.length()];
        int index = 0;
        int fp = 0;
        for(int i = 1; i< as.length()-2; i++){
            String twoCha = as.substring(i, i+2);
            if(twoCha.equals("||")){
                tmp[index] = as.substring(fp, i);
                index++;
                fp = i+2;
            }
        }
        tmp[index] = as.substring(fp,as.length());
        index++;
        String[] result = new String[index];
        for(int i = 0; i< result.length; i++){
            result[i] = tmp[i];
        }
        return result;
    }
    
    public boolean notInList(GeneInfo[] gl1, GeneInfo[] gl2, String gn1, String gn2){
		for(int i = 0; i< gl1.length; i++){
			if(gl1[i]!=null && gl1[i].name.equals(gn1) && gl2[i].name.equals(gn2)){
				return false;
			}
			if(gl1[i]!=null && gl1[i].name.equals(gn2) && gl2[i].name.equals(gn1)){
				return false;
			}
			if(gl1[i]==null){
				return true;
			}
		}
		return true;
	}

	public int getChrNumber(String chrNumber, int genomeIndex){
		int chr = -1;
		if(genomeIndex==19515){
			chr = new Integer(chrNumber.substring(18)).intValue();
		}
        
		if(genomeIndex==10997){
			chr = new Integer(chrNumber.substring(2)).intValue();
           // if(chr==0){chr=-1;}
		}
		if(genomeIndex== 9050){  // vitis
			if(chrNumber.equals("Un")){chr = 1000;}
			int l = chrNumber.length();
			if(l>7 && chrNumber.substring(l-7,l).equals("_random")){chr=1000+(new Integer(chrNumber.substring(0,l-7))).intValue();}
			if(chr==-1){
				chr = new Integer(chrNumber).intValue();
			}
			if(chr>20){chr=-1;}
		}

		
		if(genomeIndex==8400){
			chr = new Integer(chrNumber.substring(9)).intValue();
		}
		return chr;
	}
   //  int[] twoGenomeIndex = atest.getGenomeIndex(aFileName);
    
    public int[] getGenomeIndex(String aFileName){
        String[] infos = aFileName.split(Pattern.quote("."));
        
        if (infos.length < 2) {
            System.out.println("Could not get genome index for " + aFileName);
            System.exit(1);
        }
        String[] twoGenomeIndex = infos[0].split(Pattern.quote("_"));
        int[] result = new int[2];
        result[0] = Integer.parseInt(twoGenomeIndex[0]);
       
        result[1] = Integer.parseInt(twoGenomeIndex[1]);
        return result;
    }
    
	public GenePair[] removeAnEdge(GenePair[] edges){
		int index = -1;
		int weight = 101;
		for(int i = 0; i< edges.length; i++){
			if(edges[i]!=null && edges[i].weight<weight){
				index = i;
				weight = edges[i].weight;
			}
		}
		edges[index] = null;
		return edges;
	}
    
    public static String Usage() {
        return "TestGetGenomes:";
    }
}