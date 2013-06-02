opt:
	@cd f; gmake KIND=$@ OPTFLAGS=-O3
	@cd obj; gmake KIND=$@ NP=$(NP)

debug:
	@cd f; gmake KIND=$@ OPTFLAGS=-g
	@cd obj; gmake KIND=$@ NP=$(NP)

app:    
	@cd f; gmake KIND=$@ OPTFLAGS=-eA
	@cd obj; gmake KIND=$@ NP=$(NP)

test:
	@cd testf; gmake KIND=$@ OPTFLAGS=-g
	@cd obj; gmake KIND=$@ NP=$(NP)

clear:
	@rm $(HOME)/bin/higgs_$(NP)
