# RNA-seq pipeline
构建一个新的 RNA-seq pipeline  
by lianlin,zhangxuan,gaoxiaoyang

### 目录结构
```
.
├── Project #项目目录
│   ├── Analysis #工作目录
│   └── Rawdata #测试数据
└── Software #我们的 pipeline
```
### 软件下载
1. 安装 [Miniconda](https://conda.io/miniconda.html)，推荐 Python3版本
2. 设置频道
`conda config --add channels bioconda`
3. 使用Bioconda([https://bioconda.github.io/](https://bioconda.github.io/))安装 pipeline 中依赖的软件。`sh install_software.sh` or `conda install bwa`
4. 填写配置文件

### Snakemake 使用
* [官方文档](http://snakemake.readthedocs.io/en/stable/tutorial/basics.html)
* [中文简单教程](https://adoubi.life/posts/2017-11-14/)
* [例子数据](https://github.com/l0o0/snakemake-tutorial)

