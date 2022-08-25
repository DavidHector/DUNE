#ifndef DUNE_NONLINOPT_CIRCULARBUFFER_HH
#define DUNE_NONLINOPT_CIRCULARBUFFER_HH

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief A simple circular buffer
     *
     * This class provides a circular buffer that can
     * be used to store values or vectors in a sliding
     * window, with the oldest values being overwritten
     * when new entries are stored. The data type needs
     * to be default constructable, because the buffer
     * is pre-filled.
     *
     * @tparam Data type of data that should be stored
     */
    template<typename Data>
      class CircularBuffer
      {
        std::size_t offset = 0, cur_size = 0, max_size;

        Data* buffer;

        public:

        /**
         * @brief Constructor
         *
         * Allocates storage
         *
         * @param max_size_ buffer capacity
         */
        CircularBuffer(std::size_t max_size_)
          : max_size(max_size_)
        {
          buffer = new Data[max_size];
        }

        /**
         * @brief Copy constructor
         *
         * Copies storage from other circular buffer
         *
         * @param other buffer that should be copied
         */
        CircularBuffer(const CircularBuffer& other)
          : max_size(other.max_size)
        {
          buffer = new Data[max_size];
          for (std::size_t i = 0; i < max_size; i++)
            buffer[i] = other.buffer[i];
        }

        /**
         * @brief Move constructor
         *
         * Reuses storage from other circular buffer
         *
         * @param other buffer that should release storage
         */
        CircularBuffer(CircularBuffer&& other)
        : max_size(other.max_size), buffer(std::move(other.buffer))
        {
          other.buffer = nullptr;
        }

        /**
         * @brief assignment operator
         *
         * Copies storage from other circular buffer
         *
         * @param other buffer that should be copied
         *
         * @return reference to new circular buffer
         */
        CircularBuffer& operator=(const CircularBuffer& other)
        {
          max_size = other.max_size;
          delete[] buffer;
          buffer = new Data[max_size];
          for (std::size_t i = 0; i < max_size; i++)
            buffer[i] = other.buffer[i];
          return *this;
        }

        /**
         * @brief move assignement operator
         *
         * Reuses storage from other circular buffer
         *
         * @param other buffer that should release storage
         *
         * @return reference to new circular buffer
         */
        CircularBuffer& operator=(CircularBuffer&& other)
        {
          max_size = other.max_size;
          buffer = std::move(other.buffer);
          other.buffer = nullptr;
          return *this;
        }

        /**
         * @brief Destructor
         *
         * Deallocates storage
         */
        ~CircularBuffer()
        {
          delete[] buffer;
        }

        /**
         * @brief Current size
         *
         * @return number of stored objects
         */
        std::size_t size() const
        {
          return cur_size;
        }

        /**
         * @brief Check if buffer is empty
         *
         * @return true if size is zero, else false
         */
        bool empty() const
        {
          return cur_size == 0;
        }

        /**
         * @brief Rotate the circluar buffer by one
         */
        void shift()
        {
          if (cur_size < max_size)
            cur_size++;
          else
            offset = (offset + 1) % max_size;
        }

        /**
         * @brief Set size of buffer to zero
         */
        void clear()
        {
          cur_size = 0;
        }

        /**
         * @brief Reduce size of buffer by one, removing first entry
         */
        void remove_first()
        {
          if (cur_size == 0)
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::InvalidStateException,"buffer is empty");
#else
            throw std::logic_error("buffer is empty");
#endif // HAVE_DUNE_COMMON
          }

          offset = (offset + 1) % max_size;
          cur_size--;
        }

        /**
         * @brief Reduce size of buffer by one, removing last entry
         */
        void remove_last()
        {
          if (cur_size == 0)
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::InvalidStateException,"buffer is empty");
#else
            throw std::logic_error("buffer is empty");
#endif // HAVE_DUNE_COMMON
          }

          cur_size--;
        }

        /**
         * @brief Access to the i-th entry
         *
         * @return reference to requested entry
         */
        Data& operator[](int i)
        {
          if (i < 0)
            i += cur_size;

          return buffer[(offset + i) % max_size];
        }

        /**
         * @brief Read-only access to the i-th entry
         *
         * @return const reference to requested entry
         */
        const Data& operator[](int i) const
        {
          if (i < 0)
            i += cur_size;

          return buffer[(offset + i) % max_size];
        }

      };



    /**
     * @brief A circular buffer along two axes at once
     *
     * This class provides a circular buffer for storing
     * matrix entries: there are two sliding windows, one
     * for the rows and one for the columns. Shifting the
     * circular buffer shifts both the rows and columns
     * by one, forgetting the oldest row and column and
     * replacing them.
     *
     * @tparam Data type of data that should be stored
     */
    template<typename Data>
      class CircularMatrixBuffer
      {
        std::size_t offset = 0, cur_size = 0, max_size;

        Data* buffer;

        public:

        /**
         * @brief Constructor
         *
         * Allocates storage
         *
         * @param max_size_ square root of buffer capacity
         */
        CircularMatrixBuffer(std::size_t max_size_)
          : max_size(max_size_)
        {
          buffer = new Data[max_size*max_size];
        }

        /**
         * @brief Copy constructor
         *
         * Copies storage from other circular buffer
         *
         * @param other buffer that should be copied
         */
        CircularMatrixBuffer(const CircularMatrixBuffer& other)
          : max_size(other.max_size)
        {
          buffer = new Data[max_size];
          for (std::size_t i = 0; i < max_size * max_size; i++)
            buffer[i] = other.buffer[i];
        }

        /**
         * @brief Move constructor
         *
         * Reuses storage from other circular buffer
         *
         * @param other buffer that should release storage
         */
        CircularMatrixBuffer(CircularMatrixBuffer&& other)
        : max_size(other.max_size), buffer(std::move(other.buffer))
        {
          other.buffer = nullptr;
        }

        /**
         * @brief Assignment operator
         *
         * Copies storage from other circular buffer
         *
         * @param other buffer that should be copied
         *
         * @return reference to new circular buffer
         */
        CircularMatrixBuffer& operator=(const CircularMatrixBuffer& other)
        {
          max_size = other.max_size;
          delete[] buffer;
          buffer = new Data[max_size];
          for (std::size_t i = 0; i < max_size * max_size; i++)
            buffer[i] = other.buffer[i];
          return *this;
        }

        /**
         * @brief Move assignment operator
         *
         * Reuses storage from other circular buffer
         *
         * @param other buffer that should release storage
         *
         * @return reference to new circular buffer
         */
        CircularMatrixBuffer& operator=(CircularMatrixBuffer&& other)
        {
          max_size = other.max_size;
          buffer = std::move(other.buffer);
          other.buffer = nullptr;
          return *this;
        }

        /**
         * @brief Destructor
         *
         * Deallocates storage
         */
        ~CircularMatrixBuffer()
        {
          delete[] buffer;
        }

        /**
         * @brief Current size (number of rows / columns)
         *
         * @return square root of number of stored objects
         */
        std::size_t size() const
        {
          return cur_size;
        }

        /**
         * @brief Check if buffer is empty
         *
         * @return true if size is zero, else false
         */
        bool empty() const
        {
          return cur_size == 0;
        }

        /**
         * @brief Rotate the circluar buffer by one
         */
        void shift()
        {
          if (cur_size < max_size)
            cur_size++;
          else
            offset = (offset + 1) % max_size;
        }

        /**
         * @brief Set size of buffer to zero
         */
        void clear()
        {
          cur_size = 0;
        }

        /**
         * @brief Reduce size of buffer by one, removing first row / column
         */
        void remove_first()
        {
          if (cur_size == 0)
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::InvalidStateException,"buffer is empty");
#else
            throw std::logic_error("buffer is empty");
#endif // HAVE_DUNE_COMMON
          }

          offset = (offset + 1) % max_size;
          cur_size--;
        }

        /**
         * @brief Reduce size of buffer by one, removing last row / column
         */
        void remove_last()
        {
          if (cur_size == 0)
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::InvalidStateException,"buffer is empty");
#else
            throw std::logic_error("buffer is empty");
#endif // HAVE_DUNE_COMMON

          cur_size--;
        }

        /**
         * @brief Access to the (i,j)-th matrix entry
         *
         * @return reference to requested entry
         */
        Data& operator()(int i, int j)
        {
          if (i < 0)
            i += cur_size;
          if (j < 0)
            j += cur_size;

          return buffer[((offset + i) % max_size) * max_size
            + (offset + j) % max_size];
        }

        /**
         * @brief Read-only access to the (i,j)-th matrix entry
         *
         * @return const reference to requested entry
         */
        const Data& operator()(int i, int j) const
        {
          if (i < 0)
            i += cur_size;
          if (j < 0)
            j += cur_size;

          return buffer[((offset + i) % max_size) * max_size
            + (offset + j) % max_size];
        }

      };
  }
}

#endif // DUNE_NONLINOPT_CIRCULARBUFFER_HH
