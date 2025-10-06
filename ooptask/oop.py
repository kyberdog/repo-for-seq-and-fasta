class Seq:
    """
    Класс для представления биологической последовательности с заголовком.

    Атрибуты:
        header (str): Заголовок последовательности (обычно описание).
        sequence (str): Последовательность нуклеотидов или аминокислот в верхнем регистре.

    Методы:
        __str__(): Возвращает строковое представление в формате FASTA.
        __len__(): Возвращает длину последовательности.
        alphabet(): Определяет тип последовательности ("nucleotide", "protein", "likely nucleotide", 
                    "likely protein" или "unknown").
    """

    def __init__(self, header, sequence):
        self.header = header.strip()
        self.sequence = sequence.strip().upper()

    def __str__(self):
        return f">{self.header}\n{self.sequence}"

    def __len__(self):
        return len(self.sequence)

    def alphabet(self):
        """
        Определяет тип алфавита последовательности.

        Возвращает:
            str: Один из вариантов:
                - "nucleotide" — если последовательность полностью из A, C, G, T, U.
                - "protein" — если последовательность содержит только стандартные аминокислоты.
                - "likely nucleotide" — большинство символов нуклеотиды, но есть исключения.
                - "likely protein" — большинство символов аминокислоты, но есть исключения.
                - "unknown" — смешанный или неизвестный алфавит.
        """
        nucleotides = set("ACGTU")
        proteins = set("ACDEFGHIKLMNPQRSTVWY")  # стандартные однобуквенные коды аминокислот
        seqset = set(self.sequence)
        if seqset <= nucleotides:
            return "nucleotide"
        elif seqset <= proteins:
            return "protein"
        else:
            nuccount = sum(i in nucleotides for i in self.sequence)
            protcount = sum(i in proteins for i in self.sequence)
            if nuccount > protcount:
                return "likely nucleotide"
            elif protcount > nuccount:
                return "likely protein"
            else:
                return "unknown"


class FastaReader:
    """
    Класс для чтения FASTA-файлов и извлечения последовательностей.

    Атрибуты:
        filepath (str): Путь к FASTA-файлу.

    Методы:
        isfasta(): Проверяет, является ли файл FASTA по признаку первого символа '>'.
        records(): Генератор, возвращающий объекты Seq по одной записи из файла.
    """
    
    def __init__(self, filepath):
        self.filepath = filepath

    def isfasta(self):
        """
        Проверяет, начинается ли файл с символа '>', характерного для FASTA.

        Возвращает:
            bool: True, если файл существует и начинается с '>'; иначе False.
        """
        try:
            with open(self.filepath, "r") as f:
                firstchar = f.read(1)
                return firstchar == ">"
        except FileNotFoundError:
            return False

    def records(self):
        """
        Читает файл и возвращает по одной записи последовательности в виде объекта Seq.

        Возвращает:
            generator: Итератор объектов Seq для каждой записи FASTA.
        """
        header = None
        seqlines = []
        with open(self.filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield Seq(header, "".join(seqlines))
                    header = line[1:]
                    seqlines = []
                else:
                    seqlines.append(line)
            if header is not None:
                yield Seq(header, "".join(seqlines))


# Demonstration
if __name__ == "__main__":
    fastapath = "example.fasta"  # Path to FASTA file

    reader = FastaReader(fastapath)
    if not reader.isfasta():
        print("Файл не в формате FASTA или не найден")
    else:
        for record in reader.records():
            print(record)
            print(f"Length: {len(record)}")
            print(f"Alphabet: {record.alphabet()}")
            print("-" * 40)
