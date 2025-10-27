#!/bin/bash

# --- Configuração ---
# Nome do script Python que queremos testar
SCRIPT_NAME="stb-translate"
# Comando para executar o Python (pode ser "python" ou "python3")
PYTHON_CMD=""
# Diretório para guardar os ficheiros de teste
TEST_DIR="test_files"

# Cores para o output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # Sem Cor

# --- Função de Verificação ---
# Esta função verifica se o comando anterior (a conversão) foi bem-sucedido
# e se o ficheiro de saída esperado foi criado e não está vazio.
check_success() {
    # $? é o código de saída do último comando. 0 = sucesso.
    # -s "$1" verifica se o ficheiro $1 existe E não está vazio.
    if [ $? -eq 0 ] && [ -s "$1" ]; then
        echo -e " ... ${GREEN}OK${NC} (Ficheiro '$1' criado)"
    else
        echo -e " ... ${RED}FALHA${NC}"
    fi
}

# --- 1. Preparação ---
echo "--- Iniciando Testador para STB-Translate ---"

# Cria o diretório de teste (e não dá erro se já existir)
mkdir -p $TEST_DIR
echo "Diretório de teste '$TEST_DIR' preparado."

# --- 2. Criação de Ficheiros de Exemplo ---
# Vamos criar uma estrutura simples de 2 átomos de Silício (Si)
# em vários formatos de ficheiro.

echo "Gerando ficheiros de entrada de exemplo..."

# --- Exemplo: si.cif ---
cat > $TEST_DIR/si.cif << 'EOF'
data_silicon
_cell_length_a     5.43
_cell_length_b     5.43
_cell_length_c     5.43
_cell_angle_alpha  90.0
_cell_angle_beta   90.0
_cell_angle_gamma  90.0
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.00 0.00 0.00
Si2 Si 0.25 0.25 0.25
EOF

# --- Exemplo: si.poscar ---
cat > $TEST_DIR/si.poscar << 'EOF'
Silicon Unit Cell
1.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
EOF

# --- Exemplo: si.fdf (O NOVO FORMATO DE ENTRADA) ---
cat > $TEST_DIR/si.fdf << 'EOF'
SystemName          Silicon
LatticeConstant     5.43 Ang

%block LatticeVectors
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
%endblock LatticeVectors

AtomicCoordinatesFormat  Fractional

%block ChemicalSpeciesLabel
1  14  Si
%endblock ChemicalSpeciesLabel

%block AtomicCoordinatesAndAtomicSpecies
0.00 0.00 0.00 1
0.25 0.25 0.25 1
%endblock AtomicCoordinatesAndAtomicSpecies
EOF

# --- Exemplo: si.siesta (Formato antigo) ---
cat > $TEST_DIR/si.siesta << 'EOF'
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
4:
1  14  0.00 0.00 0.00
1  14  0.25 0.25 0.25
EOF

# --- Exemplo: si.xyz ---
cat > $TEST_DIR/si.xyz << 'EOF'
2
Silicon Cartesian Coords
Si 0.0000 0.0000 0.0000
Si 1.3575 1.3575 1.3575
EOF

# --- Exemplo: si.lattice (Necessário para o .xyz) ---
cat > $TEST_DIR/si.lattice << 'EOF'
5.43
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
EOF

# --- Exemplo: si.dftb ---
cat > $TEST_DIR/si.dftb << 'EOF'
2 F
Si
1 1 0.00 0.00 0.00
2 1 0.25 0.25 0.25
0.0 0.0 0.0
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
EOF

# --- Exemplo: si.xsf ---
cat > $TEST_DIR/si.xsf << 'EOF'
CRYSTAL
PRIMVEC
5.430 0.000 0.000
0.000 5.430 0.000
0.000 0.000 5.430
PRIMCOORD
2 1
14 0.0000 0.0000 0.0000
14 1.3575 1.3575 1.3575
EOF

# --- Exemplo: si.fhi ---
cat > $TEST_DIR/si.fhi << 'EOF'
lattice_vector 5.430 0.000 0.000
lattice_vector 0.000 5.430 0.000
lattice_vector 0.000 0.000 5.430
atom_frac 0.00 0.00 0.00 Si
atom_frac 0.25 0.25 0.25 Si
EOF

echo -e "${GREEN}Ficheiros de exemplo gerados em '$TEST_DIR/'${NC}"


# --- 3. Execução dos Testes ---

# Usamos --no-intro para um output mais limpo nos testes
BASE_CMD="$PYTHON_CMD $SCRIPT_NAME --no-intro"

# --- Teste de LEITORES (Todos os formatos -> CIF) ---
echo -e "\n--- Testando Leitores (Todos os formatos -> .cif) ---"

echo -n "Testando: POSCAR -> CIF"
$BASE_CMD -if poscar -i $TEST_DIR/si.poscar -of cif -o $TEST_DIR/out_from_poscar.cif
check_success "$TEST_DIR/out_from_poscar.cif"

echo -n "Testando: FDF -> CIF ${YELLOW}(NOVO!)${NC}"
$BASE_CMD -if fdf -i $TEST_DIR/si.fdf -of cif -o $TEST_DIR/out_from_fdf.cif
check_success "$TEST_DIR/out_from_fdf.cif"

echo -n "Testando: SIESTA -> CIF"
$BASE_CMD -if siesta -i $TEST_DIR/si.siesta -of cif -o $TEST_DIR/out_from_siesta.cif
check_success "$TEST_DIR/out_from_siesta.cif"

echo -n "Testando: XYZ (+lattice) -> CIF"
$BASE_CMD -if xyz -i $TEST_DIR/si.xyz --lattice $TEST_DIR/si.lattice -of cif -o $TEST_DIR/out_from_xyz.cif
check_success "$TEST_DIR/out_from_xyz.cif"

echo -n "Testando: DFTB -> CIF"
$BASE_CMD -if dftb -i $TEST_DIR/si.dftb -of cif -o $TEST_DIR/out_from_dftb.cif
check_success "$TEST_DIR/out_from_dftb.cif"

echo -n "Testando: XSF -> CIF"
$BASE_CMD -if xsf -i $TEST_DIR/si.xsf -of cif -o $TEST_DIR/out_from_xsf.cif
check_success "$TEST_DIR/out_from_xsf.cif"

echo -n "Testando: FHI -> CIF"
$BASE_CMD -if fhi -i $TEST_DIR/si.fhi -of cif -o $TEST_DIR/out_from_fhi.cif
check_success "$TEST_DIR/out_from_fhi.cif"


# --- Teste de ESCRITORES (CIF -> Todos os formatos) ---
echo -e "\n--- Testando Escritores (.cif -> Todos os formatos) ---"

echo -n "Testando: CIF -> POSCAR"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of poscar -o $TEST_DIR/out_from_cif.poscar
check_success "$TEST_DIR/out_from_cif.poscar"

echo -n "Testando: CIF -> FDF"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of fdf -o $TEST_DIR/out_from_cif.fdf
check_success "$TEST_DIR/out_from_cif.fdf"

echo -n "Testando: CIF -> XYZ"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of xyz -o $TEST_DIR/out_from_cif.xyz
check_success "$TEST_DIR/out_from_cif.xyz"

echo -n "Testando: CIF -> DFTB"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of dftb -o $TEST_DIR/out_from_cif.dftb
check_success "$TEST_DIR/out_from_cif.dftb"

echo -n "Testando: CIF -> XSF"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of xsf -o $TEST_DIR/out_from_cif.xsf
check_success "$TEST_DIR/out_from_cif.xsf"

echo -n "Testando: CIF -> FHI"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of fhi -o $TEST_DIR/out_from_cif.fhi
check_success "$TEST_DIR/out_from_cif.fhi"


# --- Teste da Flag -cf (Conversão de Coordenadas) ---
echo -e "\n--- Testando Conversão de Coordenadas (-cf) ${YELLOW}(NOVO!)${NC} ---"

# O nosso si.poscar tem coordenadas 'Direct' (Fracionárias)
echo -n "Testando: POSCAR (Direct) -> FDF (Cartesian)"
$BASE_CMD -if poscar -i $TEST_DIR/si.poscar -of fdf -o $TEST_DIR/out_cart.fdf -cf cartesian
check_success "$TEST_DIR/out_cart.fdf"

# Vamos verificar se o ficheiro de saída contém 'Ang' (indicando Cartesiano)
grep -q "Ang" $TEST_DIR/out_cart.fdf
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: Formato 'Ang' (Cartesiano) encontrado.${NC}"
else
    echo -e "   -> ${RED}Falha: Formato 'Ang' (Cartesiano) não encontrado.${NC}"
fi

# O nosso si.cif é lido como Cartesiano
echo -n "Testando: CIF (Cartesian) -> POSCAR (Direct)"
$BASE_CMD -if cif -i $TEST_DIR/si.cif -of poscar -o $TEST_DIR/out_direct.poscar -cf direct
check_success "$TEST_DIR/out_direct.poscar"

# Vamos verificar se o ficheiro de saída contém 'Direct'
grep -q "Direct" $TEST_DIR/out_direct.poscar
if [ $? -eq 0 ]; then
    echo -e "   -> ${GREEN}Verificado: Formato 'Direct' encontrado.${NC}"
else
    echo -e "   -> ${RED}Falha: Formato 'Direct' não encontrado.${NC}"
fi


# --- 4. Limpeza ---
echo -e "\n--- Testes Concluídos ---"
read -p "Deseja remover o diretório '$TEST_DIR' e todos os ficheiros de teste? (s/n) " -n 1 -r
echo # Mover para a próxima linha
if [[ $REPLY =~ ^[Ss]$ ]]; then
    echo "A limpar ficheiros de teste..."
    rm -rf $TEST_DIR
    echo -e "${GREEN}Limpeza concluída.${NC}"
else
    echo "Os ficheiros de teste foram mantidos em '$TEST_DIR/' para inspeção."
fi
