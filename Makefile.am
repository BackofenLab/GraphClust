ACLOCAL_AMFLAGS = -I m4

SUBDIRS = NSPDK_src pgma_src Array Statistics

bin_SCRIPTS = GraphClust_config.pm GraphClust.pm	\
	      ari.m \
	      matrixSignificance.R	\
	      alignCenter.pl alignCenter.sge	\
	      gc_cmsearch.pl gc_cmsearch.sge	\
	      gc_jobscript.sge \
	      gc_results_cluster.pl \
	      fasta2shrep_gspan.pl fasta2shrep_gspan.sge	\
	      fastCluster.NSPDK.sge	\
	      foldFasta.pl getAdditionalScores.pl getAdditionalScoresLocal.pl \
	      glob_results.pl glob_times.pl	\
	      graphClust_AlignRange.pl	graphClust_AlignRange.sge	\
	      graphFasta.pl gspanGroups.NSPDK.sge	\
	      MASTER_GraphClust.pl \
	      matrixClass.pl 	\
	      mloc2stockholm.pl \
	      rnaclustScores2Dist.pl scoreAln.pl	\
	      treeClass.pl clusterQuality.sh

EXTRA_DIST=   examples/config.default_global	\
	      GraphClust_config.pm GraphClust_config.pm.in GraphClust.pm	\
              ari.m \
              matrixSignificance.R	\
	      alignCenter.pl alignCenter.sge	\
	      gc_cmsearch.pl gc_cmsearch.sge	\
	      gc_jobscript.sge \
	      gc_results_cluster.pl \
	      fasta2shrep_gspan.pl fasta2shrep_gspan.sge	\
	      fastCluster.NSPDK.sge	\
	      foldFasta.pl getAdditionalScores.pl getAdditionalScoresLocal.pl \
	      glob_results.pl glob_times.pl	\
	      graphClust_AlignRange.pl	graphClust_AlignRange.sge	\
	      graphFasta.pl gspanGroups.NSPDK.sge	\
	      MASTER_GraphClust.pl	\
	      matrixClass.pl 	\
	      mloc2stockholm.pl Renviron_config	\
	      rnaclustScores2Dist.pl scoreAln.pl	\
	      treeClass.pl clusterQuality.sh

################################################################################
## following targets are for private use only

LOCAL_INSTDIR=$(HOME)/install/Software
LOCAL_INSTDIR_WEB=/scratch/rna/bisge001/Software/GraphClust
LOCAL_INSTDIR_USER=/usr/local/user
LOCAL_WEBPAGE=/$(HOME)/workspace64/website-main/web/Software/GraphClust

INSTDIR_VRNA=/usr/local/user/ViennaRNA-2.0.7

REMOTE_HOST=isis.informatik.uni-freiburg.de

.PHONY: release release-install release-install-webserver

release:
	@if ls ../.git/refs/tags | grep $(distdir) ; then \
	   echo "ERROR: tag $(distdir) exists already; use release-update to force tag or set new version in configure.ac"; \
	else \
		echo "Make release $(distdir) from current revision"; \
	  	git tag $(distdir); \
	  	echo "The new tag is committed but not pushed."; \
	  	make dist-gzip; \
	  	make dist-zip; \
	  	mv $(distdir).tar.gz Releases; \
	  	mv $(distdir).zip Releases; \
	fi
#	@if [ ! -e Releases/$(distdir).tar.gz ]; then \
#	else \
#	  echo "Release already existed."; \
#	fi

release-update:
	@echo "Make or update release $(distdir) from current revision"
	@git tag -f $(distdir)
	make dist-gzip; \
	make dist-zip; \
	mv $(distdir).tar.gz Releases; \
	mv $(distdir).zip Releases; \
	echo "The new tag is committed but not pushed.";

## generate $(distdir).tar.gz from tagged release
release-gz:
	@echo "Update to release tag $(distdir)"
	@if git checkout $(distdir) ; then \
	  echo "Generate $(distdir).tar.gz" ;\
	  make dist-gzip; make dist-zip;\
	  mv $(distdir).tar.gz Releases; \
	  mv $(distdir).zip Releases; \
	  echo "Update to tip" ;\
	  git checkout master ;\
	else \
	  false ;\
	fi


release-install: release
	cd Releases; rm -r -f $(distdir); tar xzf $(distdir).tar.gz; \
        cd $(distdir); \
	./configure --prefix=$(LOCAL_INSTDIR)/$(distdir) \
	    --enable-link-static \
            --with-vrna=$(with_vrna); \
    make; make install
    
    
release-install-user: release
	cd Releases; tar xzf $(distdir).tar.gz; \
        cd $(distdir); \
	./configure --prefix=$(LOCAL_INSTDIR_USER)/$(distdir) \
	    --enable-link-static \
            --with-vrna=$(INSTDIR_VRNA) CXX=g++44; \
    make -j8; make install
    
## compile release on remote host and install in remote host's $GLOBAL_INSTALL_DIR
release-install-remote:
	scp Releases/$(distdir).tar.gz $(REMOTE_HOST):$(LOCAL_INSTDIR_USER)/tmp
	ssh $(REMOTE_HOST) 'bash --login -c " cd $(LOCAL_INSTDIR_USER)/tmp; \
	  if test -d $(distdir) ; then rm -rf $(distdir) ; fi ;\
	  tar xzf $(distdir).tar.gz ;\
	  cd $(distdir); \
	  ./configure --prefix=$(LOCAL_INSTDIR_USER)/$(distdir) --enable-link-static --with-vrna=$(INSTDIR_VRNA);\
	  make -j8; make install ;\
	  cd .. ;\
	  cp $(distdir)/config.log $(distdir).config.log ;\
	  rm -rf $(distdir) \
	"'
	
release-install-webserver: release
	cd Releases; tar xzf $(distdir).tar.gz; \
        cd $(distdir); \
	./configure --prefix=$(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION) \
	    --enable-link-static \
            --with-vrna=$(INSTDIR_VRNA); \
	make; make install; \
	cd ..; cd ..; \
	cp ./examples/websrv_example.fa $(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION); \
	cp ./examples/websrv.config $(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION); \
	cp ./websrv_results_png.pl $(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION); \
	cp ./websrv_GraphClust.sh $(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION); \
	chmod -R o+rx $(LOCAL_INSTDIR_WEB)/$(PACKAGE_VERSION) 
	

release-website:
	cp Releases/$(distdir).tar.gz $(LOCAL_WEBPAGE)/; \
	cp Releases/$(distdir).zip $(LOCAL_WEBPAGE)/; \
	cp README $(LOCAL_WEBPAGE)/;