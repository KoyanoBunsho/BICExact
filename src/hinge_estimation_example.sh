echo 'adenylate kinase'

docker exec bel ./estimate_hinge_numbers adenylate_kinase/1ake.pdb adenylate_kinase/4ake.pdb A A bic exact
docker exec bel ./estimate_hinge_numbers adenylate_kinase/1ake.pdb adenylate_kinase/4ake.pdb A A bic lh

echo 'short linker loop'
docker exec bel ./estimate_hinge_numbers linker/1mdt.pdb linker/1ddt.pdb A A bic exact 
docker exec bel ./estimate_hinge_numbers linker/1mdt.pdb linker/1ddt.pdb A A bic lh

echo 'long linker loop'
docker exec bel ./estimate_hinge_numbers linker/1nkh.pdb linker/1fgx.pdb B A bic exact
docker exec bel ./estimate_hinge_numbers linker/1nkh.pdb linker/1fgx.pdb B A bic lh

echo 'domain swapping'
docker exec bel ./estimate_hinge_numbers 2oct.pdb pdb4n6v.pdb A A bic exact
docker exec bel ./estimate_hinge_numbers 2oct.pdb pdb4n6v.pdb A A bic lh
