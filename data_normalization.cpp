#include <getopt.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "util/func.hpp"

int main(int argc, char* argv[]) {
    //Get Option
	std::string InFileName = "./test/gene_expr.tsv";
	std::string OutFileName = "./test/gene_expr.tsv";
	std::string method = "upqt";
    
	const char * const short_opts = "hi:o:m:";
    const struct option long_opts[] =  {
    	{ "help", 0, NULL, 'h' },
		{ "input", 1, NULL, 'i' },
		{ "output", 1, NULL, 'o' },
		{ "method", 1, NULL, 'm' },
		{ NULL, 0, NULL, 0 }
    };
    int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
    while (opt != -1) {
    	switch (opt) {
    		case 'h':
                data_norm_help();    						
    			return 0;
    		case 'o':
    			OutFileName = optarg;
    			break;
    		case 'i':
    			InFileName = optarg;
    			break;
    		case 'm':
    			method = optarg;
    			break;
    		case '?':
    			data_norm_help();
    			return 0;
    		case -1:
    			break;
    		default:
    			abort();     			
		}
		opt = getopt_long( argc, argv, short_opts, long_opts, NULL );		
	}

    
    //load file
    std::vector <std::string> RowName;
	std::vector <std::vector <double> > MatrixRow;
	load(InFileName, RowName, MatrixRow);	
    //matrix inversion
    int RowNum = RowName.size();
	int ColNum = MatrixRow[0].size();
    std::vector <std::vector <double> > MatrixCol(ColNum, std::vector <double>(RowNum));
    for (int i = 0; i < RowNum; ++i) {
        for (int j = 0; j < ColNum; ++j) {
            MatrixCol[j][i] =  MatrixRow[i][j];
        }
    }
    
    //normalization
    std::vector <std::vector <double> > MatrixScaled(RowNum, std::vector <double>(ColNum));
    if (method == "upqt") {
        std::vector <double> upqts;
        for (int i = 0; i < ColNum; ++i) {
            std::vector <double> Col;
            for (int j = 0; j < RowNum; ++j) {
                if (MatrixCol[i][j] > 0) {
                    Col.push_back(MatrixCol[i][j]);
                }
            }
            std::sort(Col.begin(), Col.end());
            double upqt = sorted_vector_upperquartile(Col);
            upqts.push_back(upqt);
        }

        double upqts_mean = 0.0;
        for (int i = 0; i < upqts.size(); ++i) {
            upqts_mean = upqts_mean + upqts[i];
        }
        upqts_mean = upqts_mean / upqts.size();  
        
        std::vector <double> ratios;  
        for (int i = 0; i < ColNum; ++i) {
            ratios.push_back(upqts_mean / upqts[i]);
        }

        for (int i = 0; i < RowNum; ++i) {
            for (int j = 0; j < ColNum; ++j) {
                MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
            }
        } 
    } else if (method == "median") {
        std::vector <double> medians;
        for (int i = 0; i < ColNum; ++i) {
            std::vector <double> Col;
            for (int j = 0; j < RowNum; ++j) {
                if (MatrixCol[i][j] > 0) {
                    Col.push_back(MatrixCol[i][j]);
                }
            }
            std::sort(Col.begin(), Col.end());
            double upqt = sorted_vector_median(Col);
            medians.push_back(upqt);
        }

        double medians_mean = 0.0;
        for (int i = 0; i < medians.size(); ++i) {
            medians_mean = medians_mean + medians[i];
        }
        medians_mean = medians_mean / medians.size();  
        
        std::vector <double> ratios;  
        for (int i = 0; i < ColNum; ++i) {
            ratios.push_back(medians_mean / medians[i]);
        }
        
        for (int i = 0; i < RowNum; ++i) {
            for (int j = 0; j < ColNum; ++j) {
                MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
            }
        } 
    } else if (method == "deseq") {
        //对每一行的数据求几何平均数，如果该行的几何平均数不为零（即该行所有值都不为零），则该行每个数据除以几何平均数，得到一个比例数的矩阵。
        //比例数矩阵和原数据矩阵，列数一样。
        std::vector <std::vector <double> > MatrixRatio; 
        for (int i = 0; i < RowNum; ++i) {
            double geo_mean = geometric_mean(MatrixRow[i]);
            if (geo_mean > 0) {
                std::vector <double> item;
                for (int j = 0; j < ColNum; ++j) {
                    item.push_back((MatrixRow[i][j] / geo_mean));
                }
                MatrixRatio.push_back(item);
            } 
        }
        //转置比例数矩阵
        int MatrixRatio_RowNum = MatrixRatio.size();
        int MatrixRatio_ColNum = MatrixRatio[0].size();
        std::vector <std::vector <double> > MatrixRatio_Col(MatrixRatio_ColNum, std::vector <double>(MatrixRatio_RowNum));
        for (int i = 0; i < MatrixRatio_RowNum; ++i) {
            for (int j = 0; j < MatrixRatio_ColNum; ++j) {
                MatrixRatio_Col[j][i] =  MatrixRatio[i][j];
            }
        }
        //对比例数矩阵的每一列求中位数。
        std::vector <double> medians;
        for (int i = 0; i < MatrixRatio_ColNum; ++i) {
            std::vector <double> Col;
            for (int j = 0; j < MatrixRatio_RowNum; ++j) {
                if (MatrixRatio_Col[i][j] > 0) {
                    Col.push_back(MatrixRatio_Col[i][j]);
                }
            }
            std::sort(Col.begin(), Col.end());
            double median = sorted_vector_median(Col);
            medians.push_back(median);
        }
        //原数据矩阵的每一个值除以该列的(几何平均数的中位数)。
        double medians_mean = 0.0;
        for (int i = 0; i < medians.size(); ++i) {
            medians_mean = medians_mean + medians[i];
        }
        medians_mean = medians_mean / medians.size();  
        
        std::vector <double> ratios;  
        for (int i = 0; i < ColNum; ++i) {
            ratios.push_back(medians_mean / medians[i]);
        }

        for (int i = 0; i < RowNum; ++i) {
            for (int j = 0; j < ColNum; ++j) {
                MatrixScaled[i][j] = MatrixRow[i][j] * ratios[j];
            }
        }       
    }   
    
    //write output
    std::ofstream OutFile(OutFileName, std::ios::out);
	if (!OutFile.good()) {
		std::cerr << "Error while opening output file!" << std::endl;
		return -2;
    }
    
    for (int i = 0; i < RowNum; ++i) {
        OutFile << RowName[i];
        for (int j = 0; j < ColNum; ++j) {
            OutFile << '\t' << MatrixScaled[i][j];
        }
        OutFile << std::endl;
    }         

    return 0;
}