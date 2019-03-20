import gensim
from gensim.models import word2vec
from gensim import matutils
import logging
import tempfile



# from dna2vec git repo
class SingleKModel:
    def __init__(self, model):
        self.model = model
        self.vocab_lst = sorted(model.vocab.keys())

# from dna2vec repo
class MultiKModel:
    def __init__(self, filepath):

        # Deprecated. Use gensim.models.KeyedVectors.load_word2vec_format instead.

        self.aggregate = gensim.models.KeyedVectors.load_word2vec_format(filepath)
        #self.aggregate = word2vec.Word2Vec.load_word2vec_format(filepath, binary=False)
        self.logger = logging.Logger(self.__class__.__name__)

        vocab_lens = [len(vocab) for vocab in self.aggregate.vocab.keys()]
        self.k_low = min(vocab_lens)
        self.k_high = max(vocab_lens)
        self.vec_dim = self.aggregate.vector_size

        self.data = {}
        for k in range(self.k_low, self.k_high + 1):
            self.data[k] = self.separate_out_model(k)

    def model(self, k_len):
        """
        Use vector('ACGTA') when possible
        """
        return self.data[k_len].model

    def vector(self, vocab):
        return self.data[len(vocab)].model[vocab]

    def unitvec(self, vec):
        return matutils.unitvec(vec)

    def cosine_distance(self, vocab1, vocab2):
        return np.dot(self.unitvec(self.vector(vocab1)), self.unitvec(self.vector(vocab2)))

    def l2_norm(self, vocab):
        return np.linalg.norm(self.vector(vocab))

    def separate_out_model(self, k_len):
        vocabs = [vocab for vocab in self.aggregate.vocab.keys() if len(vocab) == k_len]
        if len(vocabs) != 4 ** k_len:
            self.logger.warn('Missing {}-mers: {} / {}'.format(k_len, len(vocabs), 4 ** k_len))

        header_str = '{} {}'.format(len(vocabs), self.vec_dim)
        with tempfile.NamedTemporaryFile(mode='w') as fptr:
            print(header_str, file=fptr)
            for vocab in vocabs:
                vec_str = ' '.join("%f" % val for val in self.aggregate[vocab])
                print('{} {}'.format(vocab, vec_str), file=fptr)
            fptr.flush()
            return SingleKModel(gensim.models.KeyedVectors.load_word2vec_format(fptr.name))
