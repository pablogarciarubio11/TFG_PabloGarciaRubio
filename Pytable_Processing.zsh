#!/bin/zsh

WORK_DIR="D:\TFG\Filtrado"
OUTPUT_DIR="D:\TFG\Pytables"
COND="Ctr"

# Crear carpeta de salida si no existe
mkdir -p "$OUTPUT_DIR"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "Error: El directorio '$WORK_DIR' no existe"
    exit 1
fi

cd "$WORK_DIR"

count=0

echo "Iniciando procesamiento de archivos PyTables..."
echo "Directorio origen: $(pwd)"
echo "Directorio destino: $OUTPUT_DIR"
echo "----------------------------------------"

for file in GSM*.h5; do
    if [[ ! -f "$file" ]]; then
        echo "No se encontraron archivos"
        continue
    fi

    GSM_NUM=""
    PREFIX=""
    output_file=""

    # Extraer GSM
    if [[ "$file" =~ (GSM[0-9]+) ]]; then
        GSM_NUM="${match[1]}"
    else
        echo "⚠️ No se pudo extraer GSM: $file"
        continue
    fi

    # Detectar tipo
    if [[ "$file" == *primir* ]]; then
        PREFIX="GlobalprimiR"
    elif [[ "$file" == *miRome* ]]; then
        PREFIX="miRome"
    else
        echo "⚠️ Tipo no reconocido: $file"
        continue
    fi

    output_file="$OUTPUT_DIR/${PREFIX}_Tabula_${GSM_NUM}_${COND}.h5"

    echo "Procesando: $file"
    echo "  → $output_file"

    rm -f "$output_file"

    if ptrepack --complevel 5 --overwrite-nodes "${file}:/matrix" "${output_file}:/matrix"; then
        echo "✅ OK"
        ((count++))
    else
        echo "❌ Error en $file"
    fi

    echo "----------------------------------------"
done

echo "Procesamiento completado: $count archivos"
echo "Archivos guardados en: $OUTPUT_DIR"

