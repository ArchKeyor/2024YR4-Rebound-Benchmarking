FROM python:3.12-slim

# Dependências de compilação para o REBOUND (núcleo em C)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Instalação de dependências
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Cópia do código-fonte e dados de teste
COPY . .

# Comando padrão: executa a suíte automatizada de testes
CMD ["pytest", "tests/", "-v"]
