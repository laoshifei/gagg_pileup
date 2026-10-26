from pathlib import Path
import os, sys
folder = Path(__file__).resolve().parent
os.environ['GAGG_NULL_PROJECT'] = str(folder.parents[2])
source = folder / 'inference-source'
sys.path.insert(0, str(source if source.exists() else folder.parents[2]))
if __name__ == '__main__':
    from pileup.analysis.alpha_pair_null_study import main
    main()
