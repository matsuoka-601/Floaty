import CopyPlugin from 'copy-webpack-plugin';
import { dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));

export default {
  mode: 'development',
  entry: './main.ts',
  output: {
    path: __dirname + '/dist',
    filename: 'bundle.js'
  },
  plugins: [
    new CopyPlugin({
      patterns: ['index.html']
    })
  ],
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true
  },
  module: {
    rules: [
      {
        test: /\.ts$/,
        use: 'ts-loader',
        exclude: /node_modules/
      },
      {
        test: /\.m?js$/,
        resolve: {
          fullySpecified: false
        }
      },
      // これ何？
      {
        test: /\.(vert|frag)$/,
        type: 'asset/source',
      }
    ]
  },
  resolve: {
    extensions: ['.ts', '.js']
  },
  devServer: {
    static: './dist',
    headers: {
      'Cross-Origin-Embedder-Policy': 'require-corp',
      'Cross-Origin-Opener-Policy': 'same-origin'
    },
    // circular dependendy を非表示
    client: {
        overlay: {
        warnings: false, 
        errors: true
        }
    } 
  }, 
};
